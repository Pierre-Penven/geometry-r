# C4 Level 3: Components of geometry_r

The crate's modules and their main dependencies. Names match the code
(`src/<module>.rs`, `src/polygon/<module>.rs`).

```mermaid
C4Component
  title Components: geometry_r
  Container_Boundary(geom, "geometry_r") {
    Component(math, "math", "Rust module", "Vec3, Mat3, wrapping, bisection")
    Component(ellipsoid, "ellipsoid", "Rust module", "WGS84, constants, ray intersection")
    Component(coords, "coords", "Rust module", "Geodetic, ECEF, LocalFrame (ENU)")
    Component(geodesic, "geodesic", "Rust module", "Great-circle, Vincenty, cross-track")
    Component(satellite, "satellite", "Rust module", "SatelliteFrame, ScanAngles")
    Component(look, "look", "Rust module", "Visibility, GSO arc, polarization skew")
    Component(footprint, "footprint", "Rust module", "Coverage contour, beam footprint, Beam")
    Component(link, "link", "Rust module", "Path loss, Doppler, PFD")
    Component(sun, "sun", "Rust module", "GMST, Sun, outage, eclipse")
    Component(polygon, "polygon", "Rust module", "Cap, Polygon, MultiPolygon, containment")
    Component(boolean, "polygon::boolean", "Rust module", "Chart / cube-face overlay")
    Component(anti, "polygon::antimeridian", "Rust module", "RFC 7946 split, LonLatBox")
    Component(planar, "polygon::planar", "Rust module", "Winding offset, orientation adapter")
    Component(geojson, "geojson", "Rust module", "Parse / write GeoJSON")
    Component(batch, "batch", "Rust module", "Batch and parallel APIs")
  }
  Container_Ext(overlay, "i_overlay", "Rust crate")
  Container_Ext(json, "serde_json", "Rust crate")

  Rel(coords, ellipsoid, "Converts coordinates with")
  Rel(satellite, coords, "Builds frames from")
  Rel(look, satellite, "Places GEO slots with")
  Rel(footprint, satellite, "Casts beam rays from")
  Rel(footprint, polygon, "Returns footprints as")
  Rel(boolean, planar, "Fixes winding and orientation with")
  Rel(anti, planar, "Fixes winding and orientation with")
  Rel(boolean, overlay, "Runs overlays through")
  Rel(anti, overlay, "Clips to the world window with")
  Rel(geojson, anti, "Writes RFC 7946 output with")
  Rel(geojson, json, "Parses JSON with")
  Rel(batch, footprint, "Tests many points with")
```

`math` and `error` are used by every module and are left out of the arrows.
