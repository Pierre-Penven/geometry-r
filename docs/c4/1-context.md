# C4 Level 1: System Context

geometry_r is a library linked into a host application. It has no runtime
dependencies outside the host process.

```mermaid
C4Context
  title System Context: geometry_r
  Person(dev, "Application developer", "Builds services that need geographic and satellite geometry")
  System_Boundary(host, "Host application") {
    System(geom, "geometry_r", "Rust library: coordinates, angles, footprints, spherical polygons, GeoJSON")
  }
  System_Ext(gis, "GeoJSON producers / consumers", "GIS tools and services exchanging polygons")
  System_Ext(eph, "Ephemeris source", "Provides satellite positions and velocities (out of scope)")
  Rel(dev, geom, "Calls the public Rust API of")
  Rel(gis, host, "Sends and receives GeoJSON polygons to/from", "GeoJSON text")
  Rel(eph, host, "Supplies satellite ECEF states to", "host-specific")
```

**Elements**

| Element | Description |
|---|---|
| Application developer | Writes the host code that calls geometry_r. |
| Host application | Any Rust program that depends on the crate. It does all I/O. |
| geometry_r | This crate. |
| GeoJSON producers / consumers | Sources and sinks of polygon data; the host passes the text or `serde_json::Value` to `geometry_r::geojson`. |
| Ephemeris source | Out of scope. The host computes satellite positions and passes them as ECEF vectors. |
