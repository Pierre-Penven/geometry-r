# C4 Level 2: Containers

geometry_r is a single container: one Rust crate compiled into the host binary.
The diagram shows it with the third-party crates it is built on.

```mermaid
C4Container
  title Containers: geometry_r inside a host process
  System_Boundary(proc, "Host process") {
    Container(host, "Host application", "Rust", "Owns I/O and calls the library")
    Container(geom, "geometry_r", "Rust crate (lib)", "Pure geometry functions")
    Container(overlay, "i_overlay", "Rust crate", "Planar polygon boolean engine")
    Container(json, "serde_json", "Rust crate", "JSON value model and parser")
    Container(rayon, "rayon", "Rust crate (optional)", "Data-parallel batch APIs, feature 'parallel'")
  }
  Rel(host, geom, "Calls", "Rust API")
  Rel(geom, overlay, "Runs planar overlays through", "Rust API")
  Rel(geom, json, "Parses and builds GeoJSON with", "serde_json::Value")
  Rel(geom, rayon, "Spreads batch work across cores with", "par_iter")
```

**Elements**

| Element | Technology | Notes |
|---|---|---|
| geometry_r | Rust 1.81 library | Modules listed in [3-components.md](3-components.md). |
| i_overlay | Rust crate 1.10 | Used only through `polygon::boolean` and `polygon::antimeridian`, behind the orientation adapter in `polygon::planar`. |
| serde_json | Rust crate 1.x | Used only by `geojson`. |
| rayon | Rust crate 1.10 (optional) | Used only by `batch::par`. |
| thiserror | Rust crate 1.x | Derives `GeoError` at compile time; not shown. |
