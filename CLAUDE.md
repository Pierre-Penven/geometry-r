# geometry-r — CLAUDE.md

A fast, general-purpose Rust library for geographic and satellite geometry:
coordinates, distances, look and scan angles, footprints and spherical polygon
operations on a round Earth.

The library has no I/O and no async code. It is pure, deterministic and
thread-safe, so any service can call it from any runtime.

Note: the library is general-purpose. Do not mention any specific application
(e.g. a satcom controller) in code, docs or crate metadata.

---

## 1. Goals and non-goals

### Goals
- **Correct on a round Earth.** Every operation must work across the antimeridian
  (±180°) and at the poles. Nothing may assume a flat lon/lat plane unless the
  function says so explicitly.
- **High performance.** Hot paths do not allocate. Types are `Copy` vectors and
  matrices. Callers can precompute frames (ENU, satellite body). Batch APIs are
  available, with optional parallel versions (`parallel` feature, rayon).
- **Explicit units.** Radians and metres inside the crate. Degrees only at the
  edges (constructors, GeoJSON). Field names carry their unit (`lat_rad`,
  `range_m`, `azimuth_rad`).
- **Interoperability.** Reads and writes GeoJSON (RFC 7946), with antimeridian
  splitting on output.

### Non-goals (for now)
- Orbit propagation beyond simple models: SGP4/TLE is a later milestone (§6).
- Rendering, tiling, or map projections for display.
- Precise astronomy (IERS EOP, nutation). Low-precision Sun and GMST models are
  enough for sun-outage and eclipse prediction.

---

## 2. Toolchain and constraints

- **Edition 2024, MSRV Rust 1.88** (`rust-version = "1.88"` in Cargo.toml),
  developed on the latest stable.
  - 1.88 is the real minimum: the `i_overlay` 4 dependency chain uses `let`
    chains (`i_shape` 1.18). CI runs the tests on 1.88. Locally, run
    `cargo +1.88.0 test --all-features` before raising a dependency.
  - The 2024 edition's resolver picks versions compatible with
    `rust-version` when a crate declares one.
- Dependencies: `i_overlay` 4.5 (planar boolean engine), `serde_json` 1,
  `thiserror` 2, and `rayon` 1.12 (optional). Dev dependency: `approx`.
- Follow the `rust-development`, `unit-testing-coverage`, `clean-code` and
  `application-documentation` skills. Before saying any work is done, run
  `make check`. It runs clippy `-D warnings`, the tests (including README
  doctests), `cargo doc`, the doc checks, and `cargo fmt --check`.
- CI: `.github/workflows/ci.yml` runs on every PR and push to `main`:
  - fmt, clippy, tests (all features and default features), rustdoc
    `-D warnings` and the doc checks, on the latest stable
  - tests on the MSRV (1.88)
  - `cargo audit`
  - `cargo llvm-cov` (lcov uploaded as an artifact)
  Keep it green before merging.
- Dependency audit: `make audit` (`cargo install cargo-audit --locked`).
  Expect a clean result.
- Coverage: `make cov` (`cargo install cargo-llvm-cov --locked` plus
  `rustup component add llvm-tools-preview`). Branch coverage needs nightly,
  so use region coverage. Keep non-test regions at 100%; the regions still
  reported are assertion failure arms in test code.
- Documentation:
  - `docs/arc42/architecture.md`, `docs/c4/*.md` (Mermaid),
    `docs/technical-design.md`, `docs/adr/`.
  - `scripts/docs/new_adr.sh "title"` creates a new ADR.
  - When you add, remove or rename a module or public type, update the C4
    component diagram and the technical-design component table.
  - README code blocks are compiled as doctests (`ReadmeDoctests` in
    `lib.rs`), so keep their printed outputs real.
- `unsafe_code = "forbid"`, `missing_docs = "warn"`.
- Errors: one `GeoError` enum (thiserror). Never call `unwrap` on input. Tests
  assert specific error variants.

### i_overlay facts (verified on 4.5)
- Points are `[f64; 2]` and shapes are `Vec<Vec<[f64; 2]>>`. Use
  `SingleFloatOverlay::overlay(&clip, OverlayRule, FillRule)`.
- In a y-up frame, **`FillRule::Positive` fills counter-clockwise contours**.
  Clockwise contours count as negative. This is the opposite of 1.x, which
  forced a y-mirroring workaround that is now gone.
- Output: outer contours are counter-clockwise, holes clockwise, and contours
  are not closed (no repeated first point). Do not rely on this silently: every
  output goes through `polygon::planar::orient_shape`.
- An empty clip set works with `OverlayRule::Subject`, which is how a self-union
  (dissolve) is done.
- Precision: floats are snapped to an adaptive i32 grid, about extent/2³⁰. This is
  about 6 mm on a cube face and about 6 cm on an 85° gnomonic chart.

---

## 3. Functional requirements

IDs are stable. Reference them in tests and commits (e.g. `GEO-COORD-03`).

### 3.1 Coordinates and frames — `coords`, `ellipsoid`, `math`
| ID | Requirement |
|----|-------------|
| GEO-COORD-01 | WGS-84 ellipsoid and custom ellipsoids/spheres; constants for mean and authalic radius, Earth rotation rate, GM, GEO radius, c. |
| GEO-COORD-02 | Geodetic → ECEF conversion. |
| GEO-COORD-03 | ECEF → geodetic, closed form (Heikkinen), no iteration. Accuracy ≤ 1 mm and ≤ 1e-11 rad from −400 m to beyond GEO; sphere special case. |
| GEO-COORD-04 | Local ENU frame (precomputed rotation): ECEF↔ENU for points and directions. |
| GEO-COORD-05 | Longitudes are normalised to (−π, π]. Latitude is validated to [−90°, 90°]. Non-finite input is rejected. |
| GEO-COORD-06 | Ray/ellipsoid intersection (optionally raised by a height) and a segment-blocked-by-Earth test. |
| GEO-COORD-07 | Rotation matrices and roll/pitch/yaw attitude composition. |

### 3.2 Distances and paths — `geodesic`
| ID | Requirement |
|----|-------------|
| GEO-DIST-01 | Great-circle central angle and distance, using the unit-vector `atan2` form. Must be accurate for antipodal points and across the antimeridian. |
| GEO-DIST-02 | Haversine distance (mean radius). |
| GEO-DIST-03 | Vincenty inverse and direct on the ellipsoid (mm accuracy). Returns `NoConvergence` for near-antipodal points; `geodesic_distance` then falls back to the sphere. *(Later: Karney's algorithm for full robustness.)* |
| GEO-DIST-04 | Initial bearing, destination point, great-circle interpolation/slerp, midpoint. All pole- and antimeridian-safe. |
| GEO-DIST-05 | Cross-track distance, along-track distance, and distance to a great-circle segment. |
| GEO-DIST-06 | Slant (3D chord) distance between ECEF points. |

### 3.3 Ground ↔ satellite angles — `look`, `satellite`
| ID | Requirement |
|----|-------------|
| GEO-ANG-01 | Look angles (azimuth from geodetic north clockwise, elevation above the ellipsoid tangent, range) from a station to an ECEF target. Single and batch versions. |
| GEO-ANG-02 | Visibility above an elevation mask. Optional per-azimuth horizon mask (terrain or buildings). |
| GEO-ANG-03 | Satellite body frame: LVLH from position and velocity, a GEO nadir frame (velocity = ω×r), explicit axes, and roll/pitch/yaw bias. Choice of geocentric or geodetic nadir. |
| GEO-ANG-04 | **Scan angles** of a ground target in the body frame, in three conventions: θ/φ (off-boresight, azimuth), u/v direction cosines, and X/Y tilt (`atan2(dx,dz)`, `atan2(dy,dz)`). Must convert both ways. |
| GEO-ANG-05 | Inverse scan: body direction → ground point (ray/ellipsoid). Returns `NoIntersection` when the ray misses. |
| GEO-ANG-06 | Off-nadir angle, Earth central angle, and sub-satellite point. Closed-form relations between nadir angle, elevation, central angle, max slant range and horizon angle. |
| GEO-ANG-07 | Topocentric angular separation between two satellites seen from a station (adjacent-satellite interference). |
| GEO-ANG-08 | Minimum separation between a pointing direction and the visible **GSO arc** (NGSO arc avoidance, ITU Art. 22 / EPFD), and the GEO longitude where it occurs. |
| GEO-ANG-09 | Polarization skew angle. General vector form for any satellite with a polarization reference vector, plus the classic GEO formula. |
| GEO-ANG-10 | Earth line-of-sight between two ECEF points with a grazing altitude (inter-satellite links, gateway ↔ LEO). |

### 3.4 Spherical polygons — `polygon`
Model: vertices are unit vectors built from (lon, lat). Edges are **great-circle
minor arcs**. Each ring bounds the region **on its left** as seen from outside
the sphere. The exterior ring is counter-clockwise and holes are clockwise.

| ID | Requirement |
|----|-------------|
| GEO-POLY-01 | Build polygons from lon/lat rings. Remove closing duplicates, consecutive duplicates and spikes. Fewer than 3 vertices gives `DegeneratePolygon`. |
| GEO-POLY-02 | Orientation policy: `SmallestArea` (default for user input: take the smaller side) or `AsGiven` (trusted, so regions larger than a hemisphere are allowed; used for computed results). |
| GEO-POLY-03 | Area in steradians or m² (authalic radius), from the sum of signed spherical triangles with a fan from the ring centroid. Perimeter. |
| GEO-POLY-04 | Point-in-polygon: bounding-cap pre-reject, then parity of arc crossings against a stored reference point known to be inside (a point just left of the longest exterior edge). Arc/arc crossings use the S2-style robust 4-orientation test. Works for any polygon size, including regions that contain a pole or antipodal points. |
| GEO-POLY-05 | Bounding cap per polygon (center and radius; full sphere for large regions). Cap union, intersects, contains. |
| GEO-POLY-06 | Boolean ops on `MultiPolygon`: **intersection, union, difference, xor**, plus `union_all` and `dissolve` (merge pieces that were split at the antimeridian). |
| GEO-POLY-07 | Antimeridian-aware lon/lat bounding box (RFC 7946 style: west > east when it crosses), including great-circle latitude bulges and pole-containing regions. |
| GEO-POLY-08 | Densify (max edge angle), circle/cap polygon builder, centroid. |
| GEO-POLY-09 | Fast paths: disjoint caps give an empty intersection, a union that is simply both inputs, and a difference equal to the subject, with no overlay call. |

**Boolean algorithm (decided):**
1. **Single gnomonic chart.** Use this when every operand is "small" (it fits in a
   cap < 90°) and the merged bounding cap is ≤ 85°. Great circles map to straight
   lines, so the planar overlay (i_overlay, `FillRule::Positive`) is exact for
   great-circle edges. The chart uses a right-handed basis `e1, e2 = c × e1` and
   preserves orientation. This covers beams, service areas and countries, so it
   is the common case.
2. **Cube-face fallback** for larger extents. Clip each ring to each of the 6
   faces with 3D Sutherland–Hodgman (4 great-circle planes per face). The
   intersection point on arc a→b is `normalize(a·|n·b| + b·|n·a|)`. Project to the
   face, where it fits in the square [−1, 1]². Then fix the **winding offset** per
   polygon and face: k = [contains(t)] − w_planar(t) for a test point t inside
   the face. Add the face square |k| times with the matching orientation. This
   handles faces fully covered by a region whose boundary never enters the face.
   Overlay per face and unproject.
3. **Dissolve seams.** If the cube result fits in one chart, self-union it there.
   Otherwise return the face pieces (correct area and membership, with visible
   seams). Seam-free stitching for results larger than a hemisphere is a roadmap
   item.

### 3.5 GeoJSON — `geojson`
| ID | Requirement |
|----|-------------|
| GEO-JSON-01 | Parse Polygon, MultiPolygon, Feature, FeatureCollection and GeometryCollection (polygonal members) into `MultiPolygon`. Reject other types with `GeoJson` errors. |
| GEO-JSON-02 | Edge interpretation option. The default is `LonLatLinear { max_step_deg }`, which follows RFC 7946 by densifying lon/lat-straight edges. The other option is `GreatCircle`. Longitude jumps over 180° are read as antimeridian crossings (shortest way). An edge exactly from −180 to +180 at constant latitude means a full parallel (RFC polar caps). |
| GEO-JSON-03 | Output that follows RFC 7946. Rings are counter-clockwise exterior and clockwise holes in lon/lat. Polygons crossing the antimeridian are split into a MultiPolygon. Rings that enclose a pole are closed through ±90°. Algorithm: insert exact great-circle crossing vertices at ±180°, unwrap longitudes, close polar rings, copy holes at −360/0/+360, clip against windows [−180+360k, 180+360k], and shift back. Polygons that need no split are written exactly, without going through the overlay. |
| GEO-JSON-04 | Optional output densification so that great-circle edges render correctly in lon/lat viewers. |
| GEO-JSON-05 | Round trip: parse(write(p)) has the same area and membership as p, within tolerance. |

### 3.6 Footprints and coverage — `footprint`
| ID | Requirement |
|----|-------------|
| GEO-FP-01 | Coverage contour of a satellite for a minimum elevation, as a polygon. Per azimuth around nadir: bisection for the limb, then bisection for elevation, on the ellipsoid. |
| GEO-FP-02 | Beam footprint polygon from body frame, boresight and beam shape (circular half-angle, or elliptical with semi-axes and orientation, defined in the tangent plane). Edge rays that miss the Earth are clamped to the limb. A boresight that misses gives `NoIntersection`. |
| GEO-FP-03 | Point-in-beam test done analytically (angular), with no polygon. Fast enough for millions of terminals. |
| GEO-FP-04 | Off-boresight angle and relative gain: Gaussian main lobe `G = Gmax − 12(θ/θ3dB)²` dB, with ITU patterns (S.672, S.1528, S.465) as an extension point. |

### 3.7 Link geometry — `link`
| ID | Requirement |
|----|-------------|
| GEO-LINK-01 | Slant range, propagation delay, free-space path loss `20·log10(4πdf/c)`. |
| GEO-LINK-02 | Range rate and Doppler shift from ECEF positions and velocities. |
| GEO-LINK-03 | Power flux density from EIRP and range: `EIRP − 10·log10(4πd²)`, for checking ITU PFD limits. |

### 3.8 Sun and time — `sun`
| ID | Requirement |
|----|-------------|
| GEO-SUN-01 | Unix time → Julian date → GMST. ECI ↔ ECEF rotation about Z (no precession or nutation; documented). |
| GEO-SUN-02 | Low-precision Sun position in ECEF (about 0.01°). |
| GEO-SUN-03 | Sun outage: Sun–satellite separation seen from a station is below a threshold. Predict outage windows over a time span. |
| GEO-SUN-04 | Satellite eclipse: cylindrical shadow, then conical umbra/penumbra (power budget). |

---

## 4. Additional operations

| Pri | Operation | Status |
|-----|-----------|--------|
| P1 | **Coverage overlap and fraction.** Area of beam ∩ service area; percentage of a polygon covered by a set of beams. | Done: `MultiPolygon::overlap_area_sr`, `MultiPolygon::coverage_fraction`. |

---

## 5. Architecture

```
src/
  lib.rs          crate docs, module list, prelude
  error.rs        GeoError, Result
  math.rs         Vec3, Mat3, angle wrapping
  ellipsoid.rs    Ellipsoid (const fn new), WGS84, constants, ray intersection
  coords.rs       Geodetic, Ecef, LocalFrame (ENU), LookAngles
  geodesic.rs     great-circle, Vincenty, cross-track
  look.rs         visibility, HorizonMask, separations, GSO arc, pol skew
  satellite.rs    SatelliteFrame, ScanAngles, nadir relations
  footprint.rs    coverage contour, beam footprints, Beam (fast point-in-beam)
  link.rs         FSPL, delay, Doppler, PFD
  sun.rs          GMST, Sun position, outage windows, eclipse
  polygon/
    mod.rs        Cap, Polygon, MultiPolygon, area, contains, cleanup, enclosing_cap
    boolean.rs    gnomonic chart + cube-face overlay, coverage fraction
    antimeridian.rs  RFC 7946 split (with collinear-vertex restoration), lon/lat bbox
  geojson.rs      parse and serialize
  batch.rs        batch, `*_into` and `par::*` (feature `parallel`) APIs
examples/         quickstart.rs, bench.rs (std::time micro-benchmarks)
tests/scenarios.rs  integration: antimeridian, poles, GEO/LEO, GeoJSON round trip, membership checks
```

Conventions:
- `Ecef` is an alias of `Vec3` (metres). `Geodetic` stores radians and metres.
- Spherical polygon vertices use **geodetic** lat/lon treated as spherical, the
  same convention as GeoJSON and S2. Footprints built from ECEF ray hits go
  through `Geodetic::from_ecef(..).to_unit_vector()` so they match.
- Azimuth runs clockwise from north in [0, 2π). Longitude is in (−π, π].
- Mark public functions with `#[must_use]` where it applies. Use `#[inline]` on
  small hot functions.

---

## 6. Milestones

1. **M1 Foundations**: done.
2. **M2 Satellite geometry**: done.
3. **M3 Spherical polygons**: done.
4. **M4 GeoJSON**: done.
5. **M5 Performance**: batch/parallel APIs and benchmark done. An allocation
   audit of the polygon paths is still open.
6. **M6 Extensions** from §4: coverage overlap/fraction done.

Open items / roadmap:
- Seam-free stitching of cube-face results larger than a hemisphere.
- Statistical benchmarks with `criterion` (now possible on the current toolchain).
- Karney geodesics (removes Vincenty's antipodal failure).
- Faster point-in-polygon for big rings (edge index or gnomonic fast path).
- ITU antenna patterns (S.465, S.672, S.1528) beyond the Gaussian main lobe.

## 6b. Implementation notes (learned the hard way)
- Upgrading i_overlay from 1.x to 4.x flipped its orientation convention.
  All orientation handling lives in `polygon::planar`, so the fix was local.
  Re-check the convention (see the i_overlay facts above) after any i_overlay
  upgrade.
- Clippy rejects `!(x > 0.0)`. Write `x.is_nan() || x <= 0.0` to keep
  rejecting NaN.
- Arc crossing (`polygon::arcs_cross`) uses a half-open rule: a vertex exactly
  on the path counts as being on the negative side. Symmetric inputs, such as
  circle vertices exactly on the equator, otherwise break crossing parity.
- Chart choice uses `enclosing_cap` (Bădoiu–Clarkson) over exterior vertices.
  Unions of per-polygon caps overestimate the extent and force the cube path
  needlessly.
- i_overlay removes vertices that are collinear in its plane. This is harmless
  in a gnomonic chart, where lines are great circles, but wrong in the lon/lat
  plane. `antimeridian::Restorer` re-inserts them.
- `rustfmt.toml` sets `max_width = 120`.

## 7. Testing requirements
- Unit tests next to the code and integration tests in `tests/`. Cover positive,
  negative and boundary cases for every public function.
- Mandatory geometric edge cases: antimeridian crossing (179.5° ↔ −179.5°), exact
  ±180°, the poles (lat ±90°, polygons enclosing a pole), antipodal points,
  degenerate rings, zero-length edges, GEO altitude, and points below the
  ellipsoid.
- Reference values: Vincenty 1975 (Flinders Peak → Buninyong), known GEO look
  angles, area of a spherical cap `2πR²(1−cos r)`, and equality of boolean-op
  areas: |A∪B| = |A| + |B| − |A∩B|.
- Property-style checks: round trips (geodetic↔ECEF, scan↔ground, GeoJSON),
  and contains consistent with intersection.
