# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project
uses [Semantic Versioning](https://semver.org/).

## [0.1.0] - 2026-11-06

First release.

### Added
- **Coordinates:**
  - WGS-84 and custom ellipsoids.
  - Geodetic ↔ ECEF (closed-form Heikkinen).
  - Local ENU frames (`LocalFrame`) and look angles.
- **Distances:** great-circle, haversine, Vincenty inverse and direct,
  bearing, destination, interpolation, cross-track, along-track and
  distance to a segment.
- **Ground → target:** visibility, horizon masks, line of sight,
  topocentric separation, separation from the geostationary arc,
  polarization skew.
- **Satellite → ground:**
  - LVLH and GEO nadir frames with attitude bias.
  - Scan angles (θ/φ, u/v, X/Y tilt) in both directions.
  - Inverse scan to the ground.
  - Nadir, elevation and central-angle relations.
- **Footprints:** coverage contour for an elevation mask, circular and
  elliptical beam footprints, analytic point-in-beam (`Beam`), Gaussian
  main-lobe gain.
- **Spherical polygons:**
  - Area, perimeter, containment, bounding cap and antimeridian-aware bbox.
  - Intersection, union, difference and xor.
  - `union_all`, `dissolve`, coverage fraction.
  - All correct across ±180° and at the poles, including regions larger
    than a hemisphere.
- **GeoJSON** (RFC 7946):
  - Reading with lon/lat-linear or great-circle edges.
  - Writing with antimeridian splitting and polar closure.
- **Link and Sun:** free-space path loss, delay, Doppler, power flux density,
  GMST, Sun position, sun-outage windows, eclipse.
- **Batch APIs**, with rayon-parallel versions behind the `parallel`
  feature.

[0.1.0]: https://github.com/Pierre-Penven/geometry-r/releases/tag/v0.1.0
