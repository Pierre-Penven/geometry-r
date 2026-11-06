//! Micro-benchmarks with `std::time` (a quick, dependency-free check against
//! the performance targets in `docs/arc42/architecture.md` §10).
//!
//! Run with `cargo run --release --example bench`.

use geometry_r::batch;
use geometry_r::coords::{ecef_to_geodetic, geodetic_to_ecef};
use geometry_r::footprint::{Beam, BeamShape, coverage_contour};
use geometry_r::polygon::{MultiPolygon, Polygon, lonlat_to_unit};
use geometry_r::prelude::*;
use geometry_r::satellite::NadirReference;
use std::hint::black_box;
use std::time::Instant;

fn bench<F: FnMut() -> R, R>(name: &str, iters: u32, mut f: F) {
    for _ in 0..iters.min(10) {
        black_box(f());
    }
    let t = Instant::now();
    for _ in 0..iters {
        black_box(f());
    }
    let per = t.elapsed().as_secs_f64() / f64::from(iters);
    let (v, unit) = if per < 1e-6 {
        (per * 1e9, "ns")
    } else if per < 1e-3 {
        (per * 1e6, "µs")
    } else {
        (per * 1e3, "ms")
    };
    println!("{name:<44} {v:>10.2} {unit}");
}

fn main() -> Result<()> {
    let g = Geodetic::from_degrees(48.85, 2.35, 35.0)?;
    let p = g.to_ecef();
    bench("geodetic_to_ecef", 5_000_000, || {
        geodetic_to_ecef(&WGS84, black_box(&g))
    });
    bench("ecef_to_geodetic (Heikkinen)", 5_000_000, || {
        ecef_to_geodetic(&WGS84, black_box(p))
    });

    let frame = LocalFrame::new(&g);
    let sat = geo_satellite_position(13f64.to_radians());
    bench("LocalFrame::look_angles", 5_000_000, || {
        frame.look_angles(black_box(sat))
    });
    let sf = SatelliteFrame::geo_nadir(sat, NadirReference::Geocentric)?;
    bench("SatelliteFrame::scan_angles", 5_000_000, || {
        sf.scan_angles(black_box(p))
    });
    let a = Geodetic::from_degrees(40.7, -74.0, 0.0)?;
    bench("vincenty_inverse", 1_000_000, || {
        vincenty_inverse(&WGS84, black_box(&g), black_box(&a))
    });

    let beam = Beam::new(Vec3::Z, &BeamShape::Circular { half_angle_rad: 0.01 })?;
    bench("Beam::contains", 10_000_000, || beam.contains(&sf, black_box(p)));
    let targets: Vec<Ecef> = (0..1_000_000)
        .map(|i| {
            Geodetic::from_degrees(
                -60.0 + 120.0 * f64::from(i % 1000) / 1000.0,
                -30.0 + 0.06 * f64::from(i / 1000),
                0.0,
            )
            .map(|g| g.to_ecef())
        })
        .collect::<Result<_>>()?;
    bench("look_angles_batch (1M targets)", 5, || {
        batch::look_angles_batch(&frame, &targets)
    });

    let cov = MultiPolygon::from(coverage_contour(sat, 10f64.to_radians(), 360)?);
    let q = lonlat_to_unit([10.0, 45.0]);
    bench("contains (360-vertex coverage)", 200_000, || cov.contains(black_box(q)));
    bench("coverage_contour (360 vertices)", 50, || {
        coverage_contour(sat, 0.17, 360)
    });

    let c1 = MultiPolygon::from(Polygon::circle(lonlat_to_unit([175.0, 10.0]), 0.2, 256)?);
    let c2 = MultiPolygon::from(Polygon::circle(lonlat_to_unit([-178.0, 12.0]), 0.2, 256)?);
    bench("intersection 256×256, antimeridian (chart)", 2_000, || {
        c1.intersection(&c2)
    });
    bench("union 256×256, antimeridian (chart)", 2_000, || c1.union(&c2));
    let g1 = MultiPolygon::from(Polygon::circle(lonlat_to_unit([0.0, 0.0]), 1.4, 720)?);
    let g2 = MultiPolygon::from(Polygon::circle(lonlat_to_unit([60.0, 0.0]), 1.4, 720)?);
    bench("union of two GEO coverages (cube faces)", 100, || g1.union(&g2));
    Ok(())
}
