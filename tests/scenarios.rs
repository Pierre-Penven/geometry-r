//! End-to-end scenarios across modules: antimeridian, poles, GEO/LEO
//! geometry, GeoJSON round trips and polygon-operation consistency.

use approx::assert_abs_diff_eq;
use geometry_r::footprint::{Beam, BeamShape, beam_footprint, coverage_contour};
use geometry_r::geojson::{self, EdgeInterpretation, ReadOptions, WriteOptions};
use geometry_r::polygon::{MultiPolygon, Polygon, lonlat_to_unit};
use geometry_r::prelude::*;
use geometry_r::satellite::NadirReference;
use std::f64::consts::PI;

const GC: ReadOptions = ReadOptions {
    edges: EdgeInterpretation::GreatCircle,
    dissolve: false,
};

/// Deterministic pseudo-random points (no external RNG dependency).
fn sample_points(n: usize, seed: u64) -> Vec<Vec3> {
    let mut s = seed;
    let mut next = || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        (s >> 11) as f64 / (1u64 << 53) as f64
    };
    (0..n)
        .map(|_| {
            let z = 2.0 * next() - 1.0;
            let t = 2.0 * PI * next();
            let r = (1.0 - z * z).sqrt();
            Vec3::new(r * t.cos(), r * t.sin(), z)
        })
        .collect()
}

/// Membership of a boolean result must agree with the operands, away from
/// boundaries.
fn check_membership(a: &MultiPolygon, b: &MultiPolygon, seed: u64) {
    let i = a.intersection(b);
    let u = a.union(b);
    let d = a.difference(b);
    let x = a.xor(b);
    for p in sample_points(400, seed) {
        let (ia, ib) = (a.contains(p), b.contains(p));
        assert_eq!(i.contains(p), ia && ib, "intersection at {p:?}");
        assert_eq!(u.contains(p), ia || ib, "union at {p:?}");
        assert_eq!(d.contains(p), ia && !ib, "difference at {p:?}");
        assert_eq!(x.contains(p), ia != ib, "xor at {p:?}");
    }
}

#[test]
fn geo_beam_coverage_of_a_service_area() {
    // A GEO satellite at 13°E with a 7×7 grid of 0.25° beams around Europe.
    let sat = geo_satellite_position(13f64.to_radians());
    let frame = SatelliteFrame::geo_nadir(sat, NadirReference::Geocentric).unwrap();
    let shape = BeamShape::Circular {
        half_angle_rad: 0.25f64.to_radians(),
    };
    let center = frame
        .direction_to(Geodetic::from_degrees(48.0, 10.0, 0.0).unwrap().to_ecef())
        .unwrap();
    let c = ScanAngles::from_direction(center);
    let mut beams = Vec::new();
    for i in -3..=3 {
        for j in -3..=3 {
            let s = ScanAngles::from_tilts(
                c.tilt_x_rad + f64::from(i) * 0.4f64.to_radians(),
                c.tilt_y_rad + f64::from(j) * 0.4f64.to_radians(),
            )
            .unwrap();
            beams.push(MultiPolygon::from(
                beam_footprint(&frame, s.direction(), &shape, 48).unwrap(),
            ));
        }
    }
    let area = MultiPolygon::from(
        Polygon::from_lonlat(&[vec![[2.0, 44.0], [18.0, 44.0], [18.0, 52.0], [2.0, 52.0]]]).unwrap(),
    );
    let f = area.coverage_fraction(&beams);
    assert!(f > 0.3 && f <= 1.0, "coverage fraction {f}");
    // Analytic beam test agrees with the footprint polygon.
    let beam = Beam::new(center, &shape).unwrap();
    let fp = beam_footprint(&frame, center, &shape, 128).unwrap();
    for p in sample_points(2000, 7) {
        let ll = geometry_r::polygon::unit_to_lonlat(p);
        if fp.cap().contains(p) {
            let ecef = Geodetic::from_degrees(ll[1], ll[0], 0.0).unwrap().to_ecef();
            let r = beam.normalized_radius(&frame, ecef).unwrap_or(f64::INFINITY);
            if (r - 1.0).abs() > 0.01 {
                assert_eq!(fp.contains(p), r < 1.0);
            }
        }
    }
}

#[test]
fn pacific_coverage_across_antimeridian_to_geojson_and_back() {
    let sat = geo_satellite_position(175f64.to_radians());
    let cov = coverage_contour(sat, 20f64.to_radians(), 360).unwrap();
    assert!(cov.contains_lonlat([-170.0, 0.0]));
    let mp = MultiPolygon::from(cov);
    let geo = geojson::to_geometry(&mp, &WriteOptions::default());
    assert_eq!(geo["type"], "MultiPolygon", "split at the antimeridian");
    // Every output longitude is within [-180, 180].
    for part in geo["coordinates"].as_array().unwrap() {
        for p in part[0].as_array().unwrap() {
            let lon = p[0].as_f64().unwrap();
            assert!((-180.0..=180.0).contains(&lon));
        }
    }
    let back = geojson::parse_value(&geo, &ReadOptions { dissolve: true, ..GC }).unwrap();
    assert_eq!(back.polygons.len(), 1);
    assert_abs_diff_eq!(back.area_sr(), mp.area_sr(), epsilon = 1e-8);
    for p in sample_points(500, 3) {
        if (back.contains(p)) != mp.contains(p) {
            // Only allowed within ~1e-6 rad of the boundary (snapping).
            let near = mp.polygons[0].exterior().iter().any(|&v| v.angle_to(p) < 0.05);
            assert!(near, "membership mismatch far from boundary at {p:?}");
        }
    }
}

#[test]
fn polar_leo_footprints_union() {
    // Three LEO satellites over the north pole region.
    let footprints: Vec<MultiPolygon> = [(0.0, 85.0), (120.0, 85.0), (-120.0, 85.0)]
        .iter()
        .map(|&(lon, lat)| {
            let sat = Geodetic::from_degrees(lat, lon, 550e3).unwrap().to_ecef();
            coverage_contour(sat, 25f64.to_radians(), 120).unwrap().into()
        })
        .collect();
    let u = MultiPolygon::union_all(&footprints);
    assert!(u.contains(Vec3::Z));
    assert!(u.area_sr() > footprints[0].area_sr());
    check_membership(&footprints[0], &footprints[1], 11);
}

#[test]
fn global_geo_coverage_uses_large_polygons() {
    // Three GEO satellites 120° apart, 5° mask: union covers most of the Earth.
    let covs: Vec<MultiPolygon> = [0.0f64, 120.0, -120.0]
        .iter()
        .map(|&lon| {
            coverage_contour(geo_satellite_position(lon.to_radians()), 5f64.to_radians(), 360)
                .unwrap()
                .into()
        })
        .collect();
    let u = MultiPolygon::union_all(&covs);
    let frac = u.area_sr() / (4.0 * PI);
    assert!(frac > 0.85 && frac < 0.97, "global fraction {frac}");
    assert!(u.contains(lonlat_to_unit([180.0, 0.0])));
    assert!(!u.contains(Vec3::Z));
    check_membership(&covs[0], &covs[1], 5);
}

#[test]
fn membership_consistency_mixed_sizes() {
    let small = MultiPolygon::from(
        Polygon::from_lonlat(&[vec![[170.0, -20.0], [-160.0, -25.0], [-165.0, 15.0], [175.0, 10.0]]]).unwrap(),
    );
    let big = MultiPolygon::from(Polygon::circle(lonlat_to_unit([-150.0, 10.0]), 100f64.to_radians(), 720).unwrap());
    check_membership(&small, &big, 99);
    let hole = MultiPolygon::from(
        Polygon::from_lonlat(&[
            vec![[-10.0, -10.0], [10.0, -10.0], [10.0, 10.0], [-10.0, 10.0]],
            vec![[-5.0, -5.0], [5.0, -5.0], [5.0, 5.0], [-5.0, 5.0]],
        ])
        .unwrap(),
    );
    let bar = MultiPolygon::from(
        Polygon::from_lonlat(&[vec![[-20.0, -2.0], [20.0, -2.0], [20.0, 2.0], [-20.0, 2.0]]]).unwrap(),
    );
    check_membership(&hole, &bar, 1234);
}

#[test]
fn geojson_feature_collection_overlap() {
    let fc = r#"{"type":"FeatureCollection","features":[
        {"type":"Feature","properties":{"id":1},"geometry":{"type":"Polygon","coordinates":[[[0,0],[10,0],[10,10],[0,10],[0,0]]]}},
        {"type":"Feature","properties":{"id":2},"geometry":{"type":"Polygon","coordinates":[[[5,5],[15,5],[15,15],[5,15],[5,5]]]}}
    ]}"#;
    let mp = geojson::parse_str(fc, &ReadOptions::default()).unwrap();
    assert_eq!(mp.polygons.len(), 2);
    let a = MultiPolygon::from(mp.polygons[0].clone());
    let b = MultiPolygon::from(mp.polygons[1].clone());
    let i = a.intersection(&b);
    assert!(i.contains_lonlat([7.5, 7.5]));
    let dissolved = mp.dissolve();
    assert_eq!(dissolved.polygons.len(), 1);
    assert_abs_diff_eq!(
        dissolved.area_sr(),
        a.area_sr() + b.area_sr() - i.area_sr(),
        epsilon = 1e-10
    );
}

#[test]
fn ground_station_to_geo_and_leo_geometry() {
    let station = Geodetic::from_degrees(-33.9, 151.2, 50.0).unwrap(); // Sydney
    let frame = LocalFrame::new(&station);
    let geo = geo_satellite_position(156f64.to_radians());
    let la = frame.look_angles(geo);
    assert!(la.elevation_rad.to_degrees() > 45.0);
    assert!(
        la.azimuth_rad.to_degrees() > 0.0 && la.azimuth_rad.to_degrees() < 30.0,
        "northern sky"
    );
    // Scan angle from the satellite back to the station, and inverse.
    let sf = SatelliteFrame::geo_nadir(geo, NadirReference::Geocentric).unwrap();
    let s = sf.scan_angles(station.to_ecef()).unwrap();
    // The station is 50 m above the ellipsoid: intersect the raised surface.
    let p = sf.ground_intersection(s.direction(), &WGS84, 50.0).unwrap();
    assert!(p.distance(station.to_ecef()) < 1e-3);
    // On the bare ellipsoid the hit is a few hundred metres away.
    let g = sf.scan_to_ground(&s).unwrap();
    assert_abs_diff_eq!(g.lat_rad, station.lat_rad, epsilon = 1e-4);
    // GSO arc separation of a GEO pointing is zero.
    let sep = geometry_r::look::gso_arc_separation(&frame, geo, 0.0).unwrap();
    assert!(sep.angle_rad < 1e-6);
    // A LEO satellite below the horizon is not visible and is blocked.
    let leo = Geodetic::from_degrees(40.0, -30.0, 550e3).unwrap().to_ecef();
    assert!(!is_visible(&frame, leo, 0.0));
    assert!(!geometry_r::look::line_of_sight(station.to_ecef(), leo, 0.0));
}
