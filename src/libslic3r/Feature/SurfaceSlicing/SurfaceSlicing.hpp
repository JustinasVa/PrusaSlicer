#ifndef libslic3r_SurfaceSlicing_hpp_
#define libslic3r_SurfaceSlicing_hpp_

#include "libslic3r/Polygon.hpp"

namespace Slic3r::Arachne {
struct ExtrusionLine;
} // namespace Slic3r::Arachne

namespace Slic3r {
class PrintRegionConfig;
} // namespace Slic3r

namespace Slic3r::Feature::SurfaceSlicing {

// Is this perimeter loop on this layer subject to the surface-slicing wall
// texture (a patterned layer of the solid/pattern cadence)?
bool should_texture(const PrintRegionConfig &config, size_t layer_idx, size_t perimeter_idx, bool is_contour);

// Displace the external perimeter polygon inward in "window" sections
// (trapezoid detours with 45-degree ramps), keeping the loop one continuous
// line. Window keep-out polygons (between the original wall line and the
// dipped line) are appended to out_infill_keepout so the caller can subtract
// them from the infill area — otherwise infill would show through the
// window recesses.
Polygon apply_surface_slicing(const Polygon &polygon, const PrintRegionConfig &config, size_t layer_idx, size_t perimeter_idx, bool is_contour, Polygons &out_infill_keepout);

// Same for an Arachne extrusion line (must be closed).
Arachne::ExtrusionLine apply_surface_slicing(const Arachne::ExtrusionLine &extrusion, const PrintRegionConfig &config, size_t layer_idx, size_t perimeter_idx, bool is_contour, Polygons &out_infill_keepout);

} // namespace Slic3r::Feature::SurfaceSlicing

#endif // libslic3r_SurfaceSlicing_hpp_
