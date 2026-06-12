#include <cmath>
#include <algorithm>
#include <vector>

#include "libslic3r/Arachne/utils/ExtrusionJunction.hpp"
#include "libslic3r/Arachne/utils/ExtrusionLine.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/PrintConfig.hpp"

#include "SurfaceSlicing.hpp"

using namespace Slic3r;

namespace Slic3r::Feature::SurfaceSlicing {

// Merge displacement samples closer than this (scaled).
static double sample_merge_eps() { return scaled<double>(0.05); }
// How far outside the original wall centerline a keep-out zone extends, to
// also catch infill that overlaps into the wall band (scaled).
static double zone_outward_margin() { return scaled<double>(0.6); }

namespace {

// The checkerboard window pattern along the arc length of one loop.
// All lengths are in scaled coordinates; u is the arc-length pattern
// coordinate (already shifted to the world-anchored origin).
struct Pattern
{
    double width;    // window length along the wall
    double spacing;  // solid wall between windows
    double depth;    // inward displacement in the window middle
    double period;   // width + spacing
    double phase;    // band phase shift (stagger)

    double ramp() const { return std::min(depth, width / 2.); }

    // Inward displacement at pattern coordinate u (trapezoid with 45-degree
    // ramps; zero in the solid sections).
    double depth_at(double u) const
    {
        double p = std::fmod(u + phase, period);
        if (p < 0.)
            p += period;
        if (p <= spacing)
            return 0.;
        return std::min(depth, std::min(p - spacing, period - p));
    }

    // Is pattern coordinate u inside a window (including its endpoints)?
    bool in_window(double u) const
    {
        double p = std::fmod(u + phase, period);
        if (p < 0.)
            p += period;
        return p >= spacing - SCALED_EPSILON;
    }

    // Append every u in (u0, u1) where the trapezoid profile changes slope.
    void breakpoints(double u0, double u1, std::vector<double> &out) const
    {
        const double r     = this->ramp();
        const double bps[] = { spacing, spacing + r, period - r, period };
        for (double b : bps) {
            // Solve (u + phase) mod period == b  for u in (u0, u1).
            double first = b - phase + std::ceil((u0 + phase - b) / period + 1e-9) * period;
            for (double u = first; u < u1 - SCALED_EPSILON; u += period)
                if (u > u0 + SCALED_EPSILON)
                    out.emplace_back(u);
        }
    }
};

Pattern make_pattern(const PrintRegionConfig &config, size_t layer_idx)
{
    Pattern p;
    p.width   = scaled<double>(config.surface_slicing_hole_width.value);
    p.spacing = scaled<double>(config.surface_slicing_hole_spacing.value);
    p.depth   = scaled<double>(config.surface_slicing_hole_depth.value);
    p.period  = p.width + p.spacing;

    const int    unit = std::max(1, config.surface_slicing_solid_layers.value + config.surface_slicing_hole_layers.value);
    const size_t band = (layer_idx - 1) / size_t(unit);
    p.phase = (band % 2) * config.surface_slicing_stagger.value * p.period;
    return p;
}

// Generic over a point sequence: produce the dipped point list and the
// keep-out zones. Points are supplied as Vec2d (scaled); the closed loop is
// pts[0..n-1] with implicit closure n-1 -> 0.
struct DipResult
{
    std::vector<Vec2d> points;       // displaced loop
    std::vector<size_t> source_segment; // index of the source segment of each point
    bool modified = false;
};

DipResult dip_closed_loop(const std::vector<Vec2d> &pts, const Pattern &pattern, const bool material_inside, Polygons &out_infill_keepout)
{
    DipResult result;
    const size_t n = pts.size();
    if (n < 3)
        return result;

    // Cumulative arc length and total loop length.
    std::vector<double> cum(n + 1, 0.);
    for (size_t i = 0; i < n; ++i)
        cum[i + 1] = cum[i] + (pts[(i + 1) % n] - pts[i]).norm();
    const double total = cum[n];
    if (total < pattern.period || total <= 0.)
        return result;

    // World-anchored pattern origin: the vertex with the largest X (ties
    // resolved by Y). Unlike the seam, this is stable layer-to-layer, so the
    // windows stack vertically without any seam-position requirement.
    size_t anchor = 0;
    for (size_t i = 1; i < n; ++i)
        if (pts[i].x() > pts[anchor].x() || (pts[i].x() == pts[anchor].x() && pts[i].y() > pts[anchor].y()))
            anchor = i;
    const double s_anchor = cum[anchor];

    // Loop orientation decides which side the enclosed region is on.
    double area2 = 0.;
    for (size_t i = 0; i < n; ++i) {
        const Vec2d &a = pts[i];
        const Vec2d &b = pts[(i + 1) % n];
        area2 += a.x() * b.y() - b.x() * a.y();
    }
    const bool ccw = area2 > 0.;

    result.points.reserve(n * 2);
    result.source_segment.reserve(n * 2);

    // Keep-out zone under construction: (outer, inner) sample pairs.
    std::vector<std::pair<Vec2d, Vec2d>> zone_run;
    auto close_zone = [&]() {
        if (zone_run.size() >= 2) {
            Polygon zone;
            zone.points.reserve(zone_run.size() * 2);
            for (const auto &pr : zone_run)
                zone.points.emplace_back(coord_t(std::llround(pr.first.x())), coord_t(std::llround(pr.first.y())));
            for (auto it = zone_run.rbegin(); it != zone_run.rend(); ++it)
                zone.points.emplace_back(coord_t(std::llround(it->second.x())), coord_t(std::llround(it->second.y())));
            zone.points.erase(std::unique(zone.points.begin(), zone.points.end()), zone.points.end());
            if (zone.points.size() >= 3) {
                if (zone.is_clockwise())
                    zone.reverse();
                out_infill_keepout.emplace_back(std::move(zone));
            }
        }
        zone_run.clear();
    };

    const double eps = sample_merge_eps();
    const double margin = zone_outward_margin();
    std::vector<double> seg_samples;
    for (size_t i = 0; i < n; ++i) {
        const Vec2d &a   = pts[i];
        const Vec2d &b   = pts[(i + 1) % n];
        const Vec2d  ab  = b - a;
        const double len = ab.norm();
        if (len <= 0.)
            continue;
        const Vec2d dir = ab / len;
        // Normal pointing into the material: for the enclosed-region side
        // use the left normal on CCW loops; flip for holes (whose material
        // lies outside the enclosed region).
        Vec2d n_mat = ccw ? Vec2d(-dir.y(), dir.x()) : Vec2d(dir.y(), -dir.x());
        if (!material_inside)
            n_mat = -n_mat;

        // Sample positions within this segment: trapezoid breakpoints plus
        // the segment's end vertex.
        const double u0 = cum[i] - s_anchor;
        const double u1 = cum[i + 1] - s_anchor;
        seg_samples.clear();
        pattern.breakpoints(u0, u1, seg_samples);
        std::sort(seg_samples.begin(), seg_samples.end());
        seg_samples.emplace_back(u1);

        for (const double u : seg_samples) {
            const double t    = (u - u0) / (u1 - u0);
            const Vec2d  base = a + ab * std::clamp(t, 0., 1.);
            const double d    = pattern.depth_at(u);
            const Vec2d  q    = base + n_mat * d;
            if (result.points.empty() || (q - result.points.back()).norm() > eps) {
                result.points.emplace_back(q);
                result.source_segment.emplace_back(i);
                if (d > SCALED_EPSILON)
                    result.modified = true;
            }
            if (pattern.in_window(u))
                zone_run.emplace_back(base - n_mat * margin, q);
            else
                close_zone();
        }
    }
    close_zone();

    if (result.points.size() < 3)
        result.modified = false;
    return result;
}

} // anonymous namespace

bool should_texture(const PrintRegionConfig &config, const size_t layer_idx, const size_t perimeter_idx, const bool is_contour)
{
    const SurfaceSlicingType type = config.surface_slicing.value;
    if (type == SurfaceSlicingType::None || layer_idx <= 0 || perimeter_idx != 0)
        return false;
    if (!is_contour && type != SurfaceSlicingType::All)
        return false;

    // Layer cadence: repeating unit of solid layers followed by patterned
    // layers, anchored at the first layer above the build plate.
    const int solid = std::max(0, config.surface_slicing_solid_layers.value);
    const int holes = std::max(1, config.surface_slicing_hole_layers.value);
    const int unit  = std::max(1, solid + holes);
    return int((layer_idx - 1) % size_t(unit)) >= solid;
}

Polygon apply_surface_slicing(const Polygon &polygon, const PrintRegionConfig &config, const size_t layer_idx, const size_t perimeter_idx, const bool is_contour, Polygons &out_infill_keepout)
{
    if (!should_texture(config, layer_idx, perimeter_idx, is_contour))
        return polygon;

    std::vector<Vec2d> pts;
    pts.reserve(polygon.points.size());
    for (const Point &pt : polygon.points)
        pts.emplace_back(pt.cast<double>());

    const Pattern   pattern = make_pattern(config, layer_idx);
    const DipResult dipped  = dip_closed_loop(pts, pattern, is_contour, out_infill_keepout);
    if (!dipped.modified)
        return polygon;

    Polygon out;
    out.points.reserve(dipped.points.size());
    for (const Vec2d &q : dipped.points)
        out.points.emplace_back(coord_t(std::llround(q.x())), coord_t(std::llround(q.y())));
    return out;
}

Arachne::ExtrusionLine apply_surface_slicing(const Arachne::ExtrusionLine &extrusion, const PrintRegionConfig &config, const size_t layer_idx, const size_t perimeter_idx, const bool is_contour, Polygons &out_infill_keepout)
{
    using namespace Slic3r::Arachne;

    if (!extrusion.is_closed || extrusion.junctions.size() < 3 ||
        !should_texture(config, layer_idx, perimeter_idx, is_contour))
        return extrusion;

    // Closed Arachne extrusions duplicate the first junction at the end.
    const bool   closed_duplicate = extrusion.junctions.front().p == extrusion.junctions.back().p;
    const size_t n                = extrusion.junctions.size() - (closed_duplicate ? 1 : 0);
    if (n < 3)
        return extrusion;

    std::vector<Vec2d> pts;
    pts.reserve(n);
    for (size_t i = 0; i < n; ++i)
        pts.emplace_back(extrusion.junctions[i].p.cast<double>());

    const Pattern   pattern = make_pattern(config, layer_idx);
    const DipResult dipped  = dip_closed_loop(pts, pattern, is_contour, out_infill_keepout);
    if (!dipped.modified)
        return extrusion;

    ExtrusionLine out(extrusion.inset_idx, extrusion.is_odd, extrusion.is_closed);
    out.junctions.reserve(dipped.points.size() + 1);
    for (size_t k = 0; k < dipped.points.size(); ++k) {
        const Vec2d &q   = dipped.points[k];
        const auto  &src = extrusion.junctions[dipped.source_segment[k]];
        out.junctions.emplace_back(Point(coord_t(std::llround(q.x())), coord_t(std::llround(q.y()))), src.w, src.perimeter_index);
    }
    if (closed_duplicate) {
        // Re-close the loop on the (possibly displaced) first point.
        out.junctions.emplace_back(out.junctions.front());
    }
    return out;
}

} // namespace Slic3r::Feature::SurfaceSlicing
