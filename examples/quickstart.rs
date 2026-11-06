//! Tour of the main APIs: look angles, scan angles, footprints, polygon
//! booleans across the antimeridian and GeoJSON output.
//!
//! Run with `cargo run --example quickstart`.

use geometry_r::footprint::{BeamShape, beam_footprint, coverage_contour};
use geometry_r::geojson::{self, WriteOptions};
use geometry_r::polygon::MultiPolygon;
use geometry_r::prelude::*;
use geometry_r::satellite::NadirReference;

fn main() -> Result<()> {
    // Ground station and a GEO satellite.
    let station = Geodetic::from_degrees(-17.7, 178.0, 10.0)?; // Fiji, near the antimeridian
    let sat = geo_satellite_position((-178.0f64).to_radians());
    let look = look_angles(&station, sat);
    println!(
        "look: az {:.2}°, el {:.2}°, range {:.0} km",
        look.azimuth_rad.to_degrees(),
        look.elevation_rad.to_degrees(),
        look.range_m / 1e3
    );

    // Scan angles of the station seen from the satellite.
    let frame = SatelliteFrame::geo_nadir(sat, NadirReference::Geocentric)?;
    let scan = frame.scan_angles(station.to_ecef())?;
    println!(
        "scan: θ {:.3}°, φ {:.1}°, u {:.5}, v {:.5}",
        scan.theta_rad.to_degrees(),
        scan.phi_rad.to_degrees(),
        scan.u,
        scan.v
    );

    // A 1° beam pointed at the station, and the 10° coverage contour.
    let beam = beam_footprint(
        &frame,
        scan.direction(),
        &BeamShape::Circular {
            half_angle_rad: 1f64.to_radians(),
        },
        64,
    )?;
    let coverage = coverage_contour(sat, 10f64.to_radians(), 180)?;
    let beam = MultiPolygon::from(beam);
    let coverage = MultiPolygon::from(coverage);
    println!(
        "beam area {:.0} km², inside coverage: {:.1} %",
        beam.area_m2() / 1e6,
        100.0 * beam.coverage_fraction(std::slice::from_ref(&coverage))
    );

    // RFC 7946 output: the beam straddles ±180° and is split automatically.
    let geometry = geojson::to_geometry(&beam, &WriteOptions::default());
    println!("beam GeoJSON type: {}", geometry["type"]);
    Ok(())
}
