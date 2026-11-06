# geometry_r

[![CI](https://github.com/Pierre-Penven/geometry-r/actions/workflows/ci.yml/badge.svg)](https://github.com/Pierre-Penven/geometry-r/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/geometry_r.svg)](https://crates.io/crates/geometry_r)
[![docs.rs](https://docs.rs/geometry_r/badge.svg)](https://docs.rs/geometry_r)

A fast Rust library for geographic and satellite geometry on a round Earth.
It is for developers who need correct results everywhere on the globe,
including across the **antimeridian** (±180°) and at the **poles**:
- coordinates and local frames
- distances and bearings
- angles between ground stations and satellites, including satellite scan
  angles
- beam and coverage footprints
- polygon intersection, union and difference
- GeoJSON input and output

| Area | What you get |
|------|--------------|
| Coordinates | WGS-84 / custom ellipsoids, geodetic ↔ ECEF (closed-form), local ENU frames |
| Distances | great-circle, haversine, Vincenty, bearing, destination, interpolation, cross-track |
| Ground → target | azimuth / elevation / range, horizon masks, line of sight, angle between two satellites, angle to the geostationary arc, polarization skew |
| Satellite → ground | LVLH and GEO body frames, attitude bias, scan angles (θ/φ, u/v, X/Y tilt) in both directions |
| Footprints | coverage area for an elevation mask, circular and elliptical beam footprints, fast point-in-beam test |
| Polygons | area, containment, bounding box, **intersection / union / difference / xor**, merge, coverage fraction |
| GeoJSON | RFC 7946 reading and writing, with automatic antimeridian splitting |
| Link and Sun | path loss, delay, Doppler, power flux density, sun-outage windows, eclipse |

## Requirements

- Linux, macOS or Windows.
- Rust **1.88** or newer, with Cargo (install from <https://rustup.rs>).
- Optional, for coverage reports: `llvm-tools-preview` and `cargo-llvm-cov`
  (see [Getting help and contributing](#getting-help-and-contributing)).

## Installation

Add the crate from crates.io:

```sh
cargo add geometry_r
```

To use multi-core batch functions, enable the `parallel` feature:

```sh
cargo add geometry_r --features parallel
```

To track unreleased changes, depend on the git repository instead:

```sh
cargo add geometry_r --git https://github.com/Pierre-Penven/geometry-r
```

API documentation: <https://docs.rs/geometry_r>.

## Configuration

The library has no configuration files or environment variables. Choices are
made per call:

| Option | Where | Default |
|---|---|---|
| `parallel` Cargo feature | `Cargo.toml` | off |
| How GeoJSON edges are read | `geojson::ReadOptions { edges, .. }` | straight lines in lon/lat (RFC 7946), densified to 1° |
| Merge pieces after reading GeoJSON | `geojson::ReadOptions { dissolve, .. }` | `false` |
| Densify edges when writing GeoJSON | `geojson::WriteOptions { densify_deg }` | `None` |
| Ring orientation of user polygons | `polygon::Orientation` | `SmallestArea` (the smaller side of each ring) |
| Nadir definition of satellite frames | `satellite::NadirReference` | `Geocentric` |

## Usage

### Look angles from a ground station to a satellite

```rust
use geometry_r::prelude::*;

let paris = Geodetic::from_degrees(48.85, 2.35, 35.0)?;
let sat = geo_satellite_position(13f64.to_radians()); // GEO at 13°E
let look = look_angles(&paris, sat);
println!(
    "az {:.1}°, el {:.1}°, range {:.0} km",
    look.azimuth_rad.to_degrees(),
    look.elevation_rad.to_degrees(),
    look.range_m / 1e3
);
// az 166.0°, el 33.0°, range 38339 km
# Ok::<(), GeoError>(())
```

### Scan angles of a ground point seen from the satellite, and back

```rust
use geometry_r::prelude::*;
use geometry_r::satellite::NadirReference;

let sat = geo_satellite_position(13f64.to_radians());
let frame = SatelliteFrame::geo_nadir(sat, NadirReference::Geocentric)?;
let target = Geodetic::from_degrees(48.85, 2.35, 0.0)?;
let scan = frame.scan_angles(target.to_ecef())?;
let back = frame.scan_to_ground(&scan)?; // same point
assert!((back.lat_deg() - 48.85).abs() < 1e-6);
# Ok::<(), GeoError>(())
```

### Polygon operations across the antimeridian

```rust
use geometry_r::prelude::*;

let a = MultiPolygon::from(Polygon::from_lonlat(&[vec![
    [170.0, -5.0], [-170.0, -5.0], [-170.0, 5.0], [170.0, 5.0],
]])?);
let b = MultiPolygon::from(Polygon::from_lonlat(&[vec![
    [175.0, -10.0], [-175.0, -10.0], [-175.0, 0.0], [175.0, 0.0],
]])?);
let both = a.intersection(&b);
assert!(both.contains_lonlat([180.0, -2.0]));
println!("overlap {:.0} km²", both.area_m2() / 1e6);
// overlap 626094 km²
# Ok::<(), GeoError>(())
```

### Beam footprint as GeoJSON

```rust
use geometry_r::prelude::*;
use geometry_r::geojson::{self, WriteOptions};
use geometry_r::satellite::NadirReference;

let sat = geo_satellite_position((-178.0f64).to_radians());
let frame = SatelliteFrame::geo_nadir(sat, NadirReference::Geocentric)?;
let fiji = Geodetic::from_degrees(-17.7, 178.0, 0.0)?;
let boresight = frame.direction_to(fiji.to_ecef())?;
let beam = beam_footprint(&frame, boresight, &BeamShape::Circular { half_angle_rad: 1f64.to_radians() }, 64)?;
let json = geojson::to_geometry(&MultiPolygon::from(beam), &WriteOptions::default());
assert_eq!(json["type"], "MultiPolygon"); // split at ±180° as RFC 7946 requires
# Ok::<(), GeoError>(())
```

### Coverage fraction

```rust
use geometry_r::prelude::*;

let area = MultiPolygon::from(Polygon::from_lonlat(&[vec![[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]]])?);
let beams = vec![
    MultiPolygon::from(Polygon::from_lonlat(&[vec![[0.0, 0.0], [5.0, 0.0], [5.0, 10.0], [0.0, 10.0]]])?),
];
let f = area.coverage_fraction(&beams);
assert!((f - 0.5).abs() < 0.01);
# Ok::<(), GeoError>(())
```

To run the bundled examples:

```sh
cargo run --example quickstart
cargo run --release --example bench
```

## Running with make

| Command | What it does |
|---|---|
| `make test` | Unit tests, integration tests, doctests (including this README) |
| `make lint` | Clippy with warnings as errors |
| `make cov` | Coverage summary with missing lines |
| `make bench` | Micro-benchmarks |
| `make docs` | API docs plus the documentation checks |
| `make check` | Lint, tests, docs and format check |

## Troubleshooting

| Problem | Cause and fix |
|---|---|
| `InvalidCoordinate("latitude ... outside [-90, 90]")` | Coordinates are swapped. GeoJSON order is `[lon, lat]`, but `Geodetic::from_degrees` takes `(lat, lon, alt)`. |
| `DegeneratePolygon` | A ring has fewer than three distinct points after duplicates and spikes are removed. |
| `NoIntersection` from footprints or `scan_to_ground` | The ray points past the Earth's limb. Check the boresight. |
| `NoConvergence` from `vincenty_inverse` | The points are nearly antipodal. Use `geodesic_distance`, which falls back to the sphere. |
| A polygon covers the "wrong side" of the Earth | It is larger than a hemisphere. Build it with `Orientation::AsGiven` and a counter-clockwise exterior. |
| A union of very large regions has straight cuts | Results larger than a hemisphere are returned as cube-face pieces. Area and containment are correct. |
| Build fails with "requires rustc 1.88 or newer" | Update Rust: `rustup update stable`. |

## Getting help and contributing

- Issues and questions: <https://github.com/Pierre-Penven/geometry-r/issues>
- Before sending a change, run `make check`. For coverage, run:

  ```sh
  rustup component add llvm-tools-preview
  cargo install cargo-llvm-cov --locked
  make cov
  ```
- Project conventions and requirements: [`CLAUDE.md`](CLAUDE.md).
- Dependency check: `cargo audit` (see the technical design, section 9).

## Licence and documentation

- Licence: MIT, see [LICENSE](LICENSE).
- Architecture: [docs/arc42/architecture.md](docs/arc42/architecture.md)
- C4 diagrams: [docs/c4/](docs/c4/)
- Technical design: [docs/technical-design.md](docs/technical-design.md)
- Decisions: [docs/adr/](docs/adr/)
