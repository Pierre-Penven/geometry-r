# geometry_r — Architecture (arc42)

This document describes geometry_r as implemented in version 0.1.0. Diagrams
live in [`../c4/`](../c4/) and decisions in [`../adr/`](../adr/). Implementation
detail is in the [technical design](../technical-design.md).

---

## 1. Introduction and Goals

geometry_r is a Rust library for geographic and satellite geometry on a round
Earth. It computes:
- coordinates and local frames
- distances and bearings
- angles between ground stations and satellites, including satellite scan angles
- beam and coverage footprints
- boolean operations on spherical polygons
- GeoJSON input and output

Every operation is correct across the antimeridian and at the poles.

**Top requirements.** The full list, with stable IDs, is in `CLAUDE.md` §3.
- `GEO-COORD-*`: geodetic ↔ ECEF and ENU conversions.
- `GEO-ANG-*`: look angles, scan angles, GSO-arc separation, polarization skew.
- `GEO-POLY-*`: polygon area, containment and booleans (intersection, union,
  difference, xor) on the sphere.
- `GEO-JSON-*`: RFC 7946 input and output with antimeridian splitting.

**Quality goals**

| Priority | Goal | Meaning |
|---|---|---|
| 1 | Geometric correctness | No flat-plane assumption: results do not change when a shape moves across ±180° or over a pole. |
| 2 | Performance | Hot paths (look angles, point in beam) take tens of nanoseconds and do not allocate. Polygon booleans take well under a millisecond for typical beams. |
| 3 | Predictability | Pure, deterministic functions. Explicit units. Typed errors (`GeoError`) instead of panics on bad input. |
| 4 | Embeddability | No I/O, no global state, no async runtime. Works from any thread or runtime. |

## 2. Constraints

| Constraint | Consequence |
|---|---|
| Rust 1.81 toolchain (`rust-version = "1.81"`) | Every dependency must build on 1.81. `geo` is excluded (see [ADR 0002](../adr/0002-planar-engine-i-overlay.md)). `criterion` is not used. No float arithmetic in `const fn`. |
| Library only | No configuration files, environment variables or deployment of its own. |
| WGS-84 / GeoJSON conventions | Coordinates are `[lon, lat]` in degrees at the boundaries and radians internally. Polygons follow RFC 7946 on output. |
| `unsafe_code = "forbid"` | Only safe Rust. |

## 3. Context and Scope

See the [C4 context diagram](../c4/1-context.md).

| Neighbour | Interface |
|---|---|
| Host application (any Rust program) | Calls the public Rust API (`geometry_r::*`), linked as a crate. |
| GeoJSON producers and consumers (GIS tools, other services) | Exchange GeoJSON text or `serde_json::Value`, which the host passes to the `geojson` module. |
| Ephemeris / orbit source | Out of scope: the host provides satellite ECEF positions and velocities. |

geometry_r does not read files, open network connections or keep state between
calls.

## 4. Solution Strategy

| Problem | Approach |
|---|---|
| Round Earth and antimeridian | Represent directions as 3D unit vectors. Longitude only appears at the boundaries ([ADR 0001](../adr/0001-unit-vector-spherical-model.md)). |
| Fast coordinate conversions | Closed-form Heikkinen ECEF → geodetic. Precomputed ENU and satellite frames: one 3×3 product per point. |
| Polygon booleans on the sphere | Gnomonic projection (great circles become lines) plus the `i_overlay` planar engine. Cube-face tiling for regions too large for one chart ([ADR 0003](../adr/0003-gnomonic-chart-and-cube-faces.md)). |
| Faces fully inside a large region | Winding-offset correction from one exact point-in-polygon test per face (`polygon::planar::push_window_offset`). |
| RFC 7946 output | Insert exact crossing vertices at ±180°, unwrap longitudes, close polar rings, clip to the world window. Re-insert vertices that i_overlay dropped as collinear in lon/lat. |
| Parallelism | Optional `parallel` feature using rayon for batch APIs. The core stays single-threaded and allocation-light. |

## 5. Building Block View

See the [C4 container](../c4/2-containers.md) and
[component](../c4/3-components.md) diagrams.

**Level 1: the crate `geometry_r`.**

**Level 2: modules.**

| Module | Responsibility |
|---|---|
| `math` | `Vec3`, `Mat3`, angle and longitude wrapping, bisection helper |
| `ellipsoid` | `Ellipsoid`, `WGS84`, physical constants, ray/ellipsoid intersection, segment blocking |
| `coords` | `Geodetic`, ECEF conversions, `LocalFrame` (ENU), `LookAngles` |
| `geodesic` | great-circle, haversine, Vincenty, bearing, destination, slerp, cross-track and along-track |
| `look` | visibility, `HorizonMask`, line of sight, topocentric separation, GSO-arc separation, polarization skew |
| `satellite` | `SatelliteFrame` (LVLH / GEO nadir / attitude), `ScanAngles`, nadir/elevation relations |
| `footprint` | coverage contour, beam footprint, `Beam` (analytic point-in-beam), Gaussian gain |
| `link` | path loss, delay, Doppler, PFD |
| `sun` | GMST, Sun position, sun-outage windows, eclipse |
| `polygon` | `Cap`, `Polygon`, `MultiPolygon`: construction, area, containment |
| `polygon::boolean` | chart and cube-face overlay, `BooleanOp`, coverage fraction |
| `polygon::antimeridian` | RFC 7946 conversion, `LonLatBox` |
| `polygon::planar` (private) | planar area, winding number, winding offset, i_overlay orientation adapter |
| `geojson` | GeoJSON parsing and writing |
| `batch` | batch and `*_into` APIs, `par::*` with feature `parallel` |

## 6. Runtime View

**Scenario A: boolean operation** (`MultiPolygon::intersection`).
1. Compare the operands' bounding caps. If they are disjoint, return the
   fast-path result.
2. `fits_one_chart`: all operands are small and their enclosing cap is ≤ 85°.
   Project every ring into one gnomonic chart, overlay with i_overlay
   (`FillRule::Positive`), orient the shapes and unproject them.
3. Otherwise, for each of the 6 cube faces:
   1. Clip the rings to the face (Sutherland–Hodgman on great-circle planes).
   2. Project them onto the face.
   3. Add the winding offset.
   4. Overlay.
4. `dissolve_seams`: merge the face pieces in one chart if they fit.

**Scenario B: GeoJSON round trip.**
1. `geojson::parse_str` reads the polygonal members.
2. Edges are densified when they are lon/lat-linear.
3. `Polygon::from_lonlat` cleans the rings and orients them.
4. To write, `geojson::to_geometry` calls `to_rfc7946` on each polygon. Polygons
   that need no split are written exactly. The others go through the
   split/clip path.

**Scenario C: beam footprint.** `beam_footprint`:
1. Build a `Beam` basis around the boresight.
2. Cast each edge ray onto WGS-84.
3. Bisect rays that miss back to the limb.
4. Build a `Polygon` from the hit points.

## 7. Deployment View

No deployment of its own: geometry_r is compiled into the host binary as a
Cargo dependency. Feature `parallel` adds a rayon thread pool inside the host
process.

## 8. Crosscutting Concepts

- **Units:** radians and metres inside the library. Degrees only in
  constructors (`from_degrees`) and GeoJSON. Field names carry the unit.
- **Error handling:** one `GeoError` enum (thiserror). Every fallible public
  function returns `Result` and documents its errors. Panics happen only when
  a proven internal invariant breaks; the `expect` message says which one.
- **Numerical robustness:**
  - `atan2` forms for angles.
  - A half-open rule for arc crossings.
  - Duplicate and spike removal in rings.
  - Snapping of output coordinates to ±180 and ±90.
- **Logging and observability:** none, since the library is pure. The host
  logs.
- **Security:** no `unsafe`, no I/O. Input is validated at the API boundary
  (non-finite values, out-of-range latitudes, malformed GeoJSON).
- **Testing:**
  - Unit tests next to the code and integration scenarios in
    `tests/scenarios.rs`.
  - Membership-consistency checks of booleans on deterministic random points.
  - Coverage measured with `cargo llvm-cov`.

## 9. Architecture Decisions

| ADR | Decision |
|---|---|
| [0001](../adr/0001-unit-vector-spherical-model.md) | Spherical model on unit vectors with great-circle edges |
| [0002](../adr/0002-planar-engine-i-overlay.md) | `i_overlay` as the planar boolean engine; no `geo` |
| [0003](../adr/0003-gnomonic-chart-and-cube-faces.md) | Gnomonic chart with cube-face fallback for spherical booleans |
| [0004](../adr/0004-pin-rust-1-81.md) | Target Rust 1.81 |

## 10. Quality Requirements

Measured baselines on the development machine (`cargo run --release --example bench`).
The owner accepted the current baselines as targets (2026-11-06): a change
that makes one of them noticeably slower is a regression.

| Scenario | Baseline | Target |
|---|---|---|
| Look angles with a precomputed frame | ~22 ns | ~22 ns |
| ECEF → geodetic | ~55 ns | ~55 ns |
| Point in beam (analytic) | ~26 ns | ~26 ns |
| 1 M look angles, one core | ~21 ms | ~21 ms |
| Intersection of two 256-vertex polygons across ±180° | ~0.4 ms | ~0.4 ms |
| Union of two GEO coverages (cube path) | ~1.7 ms | ~1.7 ms |

**Correctness scenarios** (all are tests):
- Boolean results satisfy |A∪B| = |A| + |B| − |A∩B|.
- Membership of a boolean result agrees with the operands on random points.
- Antimeridian and polar inputs give the same areas as equivalent shapes
  elsewhere.

## 11. Risks and Technical Debt

| Risk / debt | Impact | Mitigation |
|---|---|---|
| Cube-face results larger than a hemisphere keep seams along cube edges | Cosmetic splits in output; area and membership are correct | Seam-free stitching (roadmap) |
| Vincenty fails for nearly antipodal points | `NoConvergence` error | `geodesic_distance` falls back to the sphere; Karney (roadmap) |
| Planar snapping in i_overlay (about extent/2³⁰) | ~6 cm on an 85° chart, ~0.1 m on split GeoJSON output | Documented; small polygons are far more precise |
| A test point exactly on a polygon edge in the winding-offset correction | Wrong offset on one face | Test points chosen at "irrational" coordinates; very unlikely |
| Point-in-polygon is O(n) per query | Slow for very large rings | Edge index (roadmap) |
| Rust 1.81 pin | Older dependency versions | Revisit when the toolchain is upgraded |
| RUSTSEC-2025-0165: `i_tree` 0.8.3 (via `i_overlay` 1.10) exposes unsound safe accessors (`Tree::node`, `Tree::mut_node`) | Undefined behaviour only if called with an out-of-range index; geometry_r never calls `i_tree`, and `i_overlay` passes internal indices | Accepted while on Rust 1.81: the fix (`i_tree` ≥ 0.10) needs `i_overlay` ≥ 2, whose dependencies need Rust ≥ 1.85 |

## 12. Glossary

| Term | Definition |
|---|---|
| ECEF | Earth-Centred Earth-Fixed Cartesian frame (metres). |
| ENU | Local East-North-Up frame at a ground point. |
| LVLH | Local-Vertical Local-Horizontal satellite frame: Z to nadir, X roughly along the velocity. |
| Scan angles | Direction of a target in a satellite body frame: θ/φ, u/v direction cosines, or X/Y tilts. |
| Nadir angle | Angle at the satellite between nadir and a target. |
| Central angle | Angle at the Earth's centre between two points. |
| GSO arc | The geostationary orbit, seen as an arc in a station's sky. |
| Antimeridian | The ±180° meridian. |
| Gnomonic projection | Central projection onto a tangent plane; maps great circles to straight lines. |
| Cap | Spherical cap: all points within an angle of a centre direction. |
| RFC 7946 | The GeoJSON specification. |
| Winding offset | Constant correction of the planar fill count inside a clipped window. |
