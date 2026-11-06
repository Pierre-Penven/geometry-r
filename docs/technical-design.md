# geometry_r — Technical Design

Implementation detail for engineers changing or embedding the library. For the
big picture, read the [arc42 document](arc42/architecture.md) and the
[C4 diagrams](c4/).

## 1. Overview

geometry_r is a pure Rust library. Main design choices:
- **Unit-vector spherical model.** Polygons and great-circle functions work on
  3D unit vectors, so there are no antimeridian or pole special cases
  ([ADR 0001](adr/0001-unit-vector-spherical-model.md)).
- **Ellipsoid where it matters.** ECEF conversions, look angles, ray casting
  and footprints use WGS-84.
- **Planar engine through charts.** Polygon booleans project to a gnomonic
  chart (or to cube faces) and run `i_overlay`
  ([ADR 0002](adr/0002-planar-engine-i-overlay.md),
  [ADR 0003](adr/0003-gnomonic-chart-and-cube-faces.md)).
- **Precompute once, apply many times.** `LocalFrame`, `SatelliteFrame` and
  `Beam` hold rotation matrices, so per-point work is a few dot products.

## 2. Components

| Module | Main items | Depends on |
|---|---|---|
| `math` | `Vec3`, `Mat3`, `wrap_pi`, `wrap_two_pi`, `wrap_lon_deg(_closed)`, `wrap_lon_delta_deg`, `bisect_last_true`, `bisect_steps_for` | — |
| `error` | `GeoError`, `Result` | thiserror |
| `ellipsoid` | `Ellipsoid` (`new`, `sphere`, radii, `normal`, `ray_intersection`, `segment_blocked`), `WGS84`, constants (`EARTH_*`, `GEO_RADIUS_M`, `SPEED_OF_LIGHT_M_S`) | math |
| `coords` | `Geodetic`, `Ecef`, `geodetic_to_ecef`, `ecef_to_geodetic` (Heikkinen), `LocalFrame`, `LookAngles`, `enu_to_look` | ellipsoid |
| `geodesic` | `central_angle`, `great_circle_distance`, `haversine_distance`, `initial_bearing`, `destination`, `interpolate`, `midpoint`, `slerp`, `cross_track_distance`, `along_track_distance`, `distance_to_segment`, `vincenty_inverse`, `vincenty_direct`, `geodesic_distance` | coords |
| `satellite` | `SatelliteFrame` (`lvlh`, `geo_nadir`, `from_axes`, `with_attitude`, `scan_angles`, `ground_intersection`, `scan_to_ground`), `ScanAngles`, `NadirReference`, nadir/elevation relations, `geo_satellite_position`, `sub_satellite_point` | coords |
| `look` | `look_angles`, `is_visible`, `HorizonMask`, `line_of_sight`, `topocentric_separation`, `gso_arc_separation`, `polarization_skew`, `geo_polarization_skew` | satellite |
| `footprint` | `BeamShape`, `Beam`, `point_in_beam`, `beam_normalized_radius`, `off_boresight_angle`, `beam_footprint`, `coverage_contour`, `gaussian_gain_db` | satellite, polygon |
| `link` | `slant_range`, `propagation_delay_s`, `free_space_path_loss_db`, `power_flux_density_dbw_m2`, `range_rate`, `doppler_shift_hz` | math |
| `sun` | `unix_to_julian_date`, `gmst_rad`, `eci_to_ecef(_matrix)`, `ecef_to_eci`, `sun_position_eci/_ecef`, `sun_satellite_separation`, `is_sun_outage`, `sun_outage_windows`, `EclipseState`, `eclipse_state`, `in_cylindrical_shadow` | math |
| `polygon` | `LonLat`, `Cap`, `Orientation`, `Polygon`, `MultiPolygon`, `lonlat_to_unit`, `unit_to_lonlat`, `ring_signed_area` | geodesic |
| `polygon::boolean` | `BooleanOp`; `MultiPolygon::{boolean, intersection, union, difference, xor, dissolve, union_all, overlap_area_sr, coverage_fraction}` | planar, i_overlay |
| `polygon::antimeridian` | `to_rfc7946`, `LonLatPolygon`, `LonLatBox`, `bounding_box` | planar, i_overlay |
| `polygon::planar` (private) | `signed_area`, `winding_number`, `push_window_offset`, `orient_shape` | — |
| `geojson` | `EdgeInterpretation`, `ReadOptions`, `WriteOptions`, `parse_str`, `parse_value`, `to_geometry`, `to_feature`, `to_string`, `to_lonlat_polygons` | antimeridian, serde_json |
| `batch` | `geodetic_to_ecef_batch`, `ecef_to_geodetic_batch`, `look_angles_batch`, `look_angles_into`, `scan_angles_batch`, `points_in_beam`, `contains_batch`; `par::*` | all above, rayon (opt.) |

## 3. Data Model

There is no persistence. The core value types:

| Type | Fields | Invariants |
|---|---|---|
| `Vec3` (= `Ecef`) | `x, y, z: f64` | ECEF positions are in metres; directions are unit vectors where documented |
| `Geodetic` | `lat_rad`, `lon_rad`, `alt_m` | lat ∈ [−π/2, π/2]; lon ∈ (−π, π]; all finite (checked by `from_degrees` / `from_radians`) |
| `LookAngles` | `azimuth_rad` [0, 2π) from north, clockwise; `elevation_rad`; `range_m` | — |
| `LocalFrame` | `origin`, `geodetic`, `rot` (rows E, N, U) | `rot` is orthonormal |
| `SatelliteFrame` | `position`, `rot` (rows body X, Y, Z in ECEF) | `rot` is orthonormal; Z is the boresight / nadir |
| `ScanAngles` | `theta_rad`, `phi_rad`, `u`, `v`, `tilt_x_rad`, `tilt_y_rad` | all describe the same unit direction |
| `Cap` | `center` (unit), `radius_rad` ∈ [0, π] | — |
| `Polygon` | `rings: Vec<Vec<Vec3>>` (private), `cap`, `small`, `inside_ref` | ring 0 is the exterior (CCW from outside), the others are holes (CW); rings are open with ≥ 3 distinct vertices; `inside_ref` lies inside the polygon |
| `MultiPolygon` | `polygons: Vec<Polygon>` | members assumed non-overlapping (use `dissolve` otherwise) |
| `LonLatBox` | `west, south, east, north` (degrees) | `west > east` when crossing ±180° |

## 4. Interfaces

The public Rust API, documented with rustdoc (`cargo doc --open`). Typical
calls:

```rust
use geometry_r::prelude::*;

// Look angles from a station to a GEO satellite.
let station = Geodetic::from_degrees(48.85, 2.35, 35.0)?;
let la = look_angles(&station, geo_satellite_position(13f64.to_radians()));

// Scan angles and the inverse.
let frame = SatelliteFrame::geo_nadir(geo_satellite_position(0.0), geometry_r::satellite::NadirReference::Geocentric)?;
let s = frame.scan_angles(station.to_ecef())?;
let ground = frame.scan_to_ground(&s)?;

// Polygons and GeoJSON.
let a = MultiPolygon::from(Polygon::from_lonlat(&[vec![[170.0, -5.0], [-170.0, -5.0], [-170.0, 5.0], [170.0, 5.0]]])?);
let json = geometry_r::geojson::to_string(&a, &Default::default()); // MultiPolygon split at ±180°
# Ok::<(), GeoError>(())
```

**GeoJSON input:**
- Accepted: `Polygon`, `MultiPolygon`, `Feature` (a null geometry yields no
  polygon), `FeatureCollection`, and `GeometryCollection` with polygonal
  members.
- Positions may carry a third value (altitude), which is ignored.
- Any other type returns `GeoError::GeoJson`.

**GeoJSON output:**
- A `Polygon` when there is exactly one part, otherwise a `MultiPolygon`.
- Closed rings, exterior counter-clockwise.

## 5. Key Flows

### Boolean operation

`MultiPolygon::boolean` → `overlay_sets`:

1. **Fast path.** If the bounding caps are disjoint, return immediately:
   - intersection: empty
   - difference: the subject
   - union and xor: both inputs together
2. **Chart path.** If `fits_one_chart` succeeds:
   1. `Chart::project_polygons` maps every ring to mirrored-y gnomonic
      coordinates.
   2. i_overlay runs with `FillRule::Positive`.
   3. `shape_to_polygon` calls `orient_shape`, unprojects, runs `clean_ring`,
      and drops collapsed shapes.
3. **Cube path** (`overlay_on_cube_faces`). For each face whose cap meets an
   operand:
   1. `clip_ring_to_face`: Sutherland–Hodgman against four planes.
   2. `push_window_offset`: fix the fill count with `poly.contains(test
      point)`.
   3. Overlay.

   Then `dissolve_seams` merges the pieces if they fit in one chart.

**Error paths:** none. Booleans are infallible. Degenerate output pieces are
dropped.

### RFC 7946 output

`to_rfc7946`:
1. Optionally densify the polygon.
2. `unwrap_ring` each ring: insert crossings at ±180°, unwrap, close polar
   rings.
3. **Exact path.** Small polygons that neither cross ±180° nor turn around a
   pole are written vertex for vertex.
4. **Split path** (`split_with_overlay`):
   1. `shifted_copies` makes ±360° copies.
   2. Apply the winding offset in the world window.
   3. Intersect with the world rectangle.
   4. `orient_shape`.
   5. `Restorer::restore` re-inserts collinear vertices.
   6. Unmirror and snap to ±180 / ±90.

### Coverage contour

For each azimuth around nadir:
1. Bisect for the limb (ray hits or misses the ellipsoid).
2. Bisect for the nadir angle whose ground elevation equals the mask.
3. Convert the hit point to a unit vector.

## 6. Configuration

There are no runtime settings. Compile-time options:

| Cargo feature | Default | Effect |
|---|---|---|
| `parallel` | off | Enables `batch::par::*` (rayon) |

Per-call options are plain structs: `geojson::ReadOptions`, `WriteOptions` and
`EdgeInterpretation`; `polygon::Orientation`; `satellite::NadirReference`.

## 7. Build, Test, and Deploy

| Task | Command |
|---|---|
| Build | `cargo build --release` |
| Format / lint | `make fmt`, `make lint` (clippy `-D warnings`) |
| Tests | `make test` (unit + `tests/scenarios.rs` + doctests) |
| Coverage | `make cov` (needs `llvm-tools-preview` and `cargo-llvm-cov 0.6.21`) |
| Benchmarks | `make bench` |
| API docs | `make docs` (also runs the doc checks) |
| Everything | `make check` |
| CI | GitHub Actions (`.github/workflows/ci.yml`): every PR and push to `main` runs the 1.81 checks, tests on latest stable, `cargo audit` and coverage |

**Deploy:** the host adds the crate as a Cargo dependency, from a local path
or from git (`https://github.com/Pierre-Penven/geometry-r`). It is not
published to a registry. TODO: verify the long-term distribution channel.

## 8. Operations

The library produces no logs or metrics, and has no health checks; the host
owns those. Failure modes surface as `GeoError` values:

| Error | Typical cause | Fix |
|---|---|---|
| `InvalidCoordinate` | latitude outside ±90°, NaN or ∞ | Validate or clean the input |
| `InvalidArgument` | bad parameter (n < 3, beam angle out of range, empty mask, degenerate frame) | Fix the call |
| `DegeneratePolygon` | a ring with fewer than 3 distinct vertices after cleaning | Drop or repair the ring |
| `NoConvergence` | Vincenty inverse for nearly antipodal points | Use `geodesic_distance` (sphere fallback) |
| `NoIntersection` | the boresight or scan ray misses the Earth | Check pointing; expected past the limb |
| `GeoJson` | malformed JSON or a non-polygonal geometry | Fix the input |

## 9. Security

- `unsafe_code = "forbid"`. No I/O, no network, no file system access.
- Inputs are validated at the boundary: finite numbers, latitude range,
  GeoJSON structure.
- Resource use grows linearly with input size (vertex counts). The library
  sets no hard limits. Hosts that accept untrusted GeoJSON should cap the
  payload size and vertex count before parsing.
- Dependencies: `i_overlay`, `serde_json`, `thiserror`, and optionally
  `rayon`. Run `cargo audit` to check them against the RustSec database.
  - The last audit (2026-11-06) found 0 vulnerabilities and 1 informational
    warning: RUSTSEC-2025-0165, `i_tree` 0.8.3 "unsound", pulled in by
    `i_overlay` 1.10.
  - geometry_r never calls the affected functions. The fix needs a newer
    toolchain (see the risks in the arc42 document).
  - The tool is `cargo-audit` 0.22.2, built with a side toolchain
    (`cargo +1.88.0 install cargo-audit --locked`) because the version that
    builds on 1.81 (0.21.2) cannot parse current advisories that use CVSS 4.0
    scores.

## 10. Known Limitations

- Results larger than a hemisphere from the cube path keep seams along cube
  edges. Area and membership are correct.
- The polygon model is spherical (geodetic latitude treated as spherical).
  Areas use the authalic radius.
- Planar snapping:
  - about 6 cm on an 85° chart
  - about 10 cm on split GeoJSON output
- Vincenty fails for nearly antipodal points.
- Sun and Earth rotation are low precision (about 0.01°; no precession,
  nutation or UT1).
- Point-in-polygon is O(n) per query: there is no spatial index.
- Only a Gaussian main-lobe antenna pattern is provided.
