//! Batch APIs: one call for many points, reusing precomputed frames.
//!
//! The `*_into` variants write into caller-provided buffers and never
//! allocate. With the `parallel` feature, `par_*` variants spread the work
//! over all cores with rayon.

use crate::coords::{ecef_to_geodetic, geodetic_to_ecef, Ecef, Geodetic, LocalFrame, LookAngles};
use crate::ellipsoid::Ellipsoid;
use crate::error::{GeoError, Result};
use crate::footprint::{Beam, BeamShape};
use crate::math::Vec3;
use crate::polygon::MultiPolygon;
use crate::satellite::{SatelliteFrame, ScanAngles};

fn check_len(a: usize, b: usize) -> Result<()> {
    if a == b {
        Ok(())
    } else {
        Err(GeoError::InvalidArgument(format!(
            "output length {b} != input length {a}"
        )))
    }
}

/// Geodetic → ECEF for many points.
#[must_use]
pub fn geodetic_to_ecef_batch(ell: &Ellipsoid, points: &[Geodetic]) -> Vec<Ecef> {
    points.iter().map(|g| geodetic_to_ecef(ell, g)).collect()
}

/// ECEF → geodetic for many points.
#[must_use]
pub fn ecef_to_geodetic_batch(ell: &Ellipsoid, points: &[Ecef]) -> Vec<Geodetic> {
    points.iter().map(|&p| ecef_to_geodetic(ell, p)).collect()
}

/// Look angles from one station to many targets.
#[must_use]
pub fn look_angles_batch(frame: &LocalFrame, targets: &[Ecef]) -> Vec<LookAngles> {
    targets.iter().map(|&t| frame.look_angles(t)).collect()
}

/// Look angles from one station to many targets, into a buffer.
///
/// # Errors
/// [`GeoError::InvalidArgument`] when the buffer length differs.
pub fn look_angles_into(frame: &LocalFrame, targets: &[Ecef], out: &mut [LookAngles]) -> Result<()> {
    check_len(targets.len(), out.len())?;
    for (o, &t) in out.iter_mut().zip(targets) {
        *o = frame.look_angles(t);
    }
    Ok(())
}

/// Scan angles of many targets from one platform (`None` for a target at
/// the platform position).
#[must_use]
pub fn scan_angles_batch(frame: &SatelliteFrame, targets: &[Ecef]) -> Vec<Option<ScanAngles>> {
    targets.iter().map(|&t| frame.scan_angles(t).ok()).collect()
}

/// Analytic point-in-beam test for many targets.
///
/// # Errors
/// [`GeoError::InvalidArgument`] for an invalid beam.
pub fn points_in_beam(
    frame: &SatelliteFrame,
    boresight_body: Vec3,
    shape: &BeamShape,
    targets: &[Ecef],
) -> Result<Vec<bool>> {
    let beam = Beam::new(boresight_body, shape)?;
    Ok(targets.iter().map(|&t| beam.contains(frame, t)).collect())
}

/// Point-in-multipolygon for many unit vectors.
#[must_use]
pub fn contains_batch(mp: &MultiPolygon, points: &[Vec3]) -> Vec<bool> {
    points.iter().map(|&p| mp.contains(p)).collect()
}

/// Parallel versions (feature `parallel`).
#[cfg(feature = "parallel")]
pub mod par {
    use super::*;
    use rayon::prelude::*;

    /// Parallel [`super::look_angles_batch`].
    #[must_use]
    pub fn look_angles_batch(frame: &LocalFrame, targets: &[Ecef]) -> Vec<LookAngles> {
        targets.par_iter().map(|&t| frame.look_angles(t)).collect()
    }

    /// Parallel [`super::geodetic_to_ecef_batch`].
    #[must_use]
    pub fn geodetic_to_ecef_batch(ell: &Ellipsoid, points: &[Geodetic]) -> Vec<Ecef> {
        points.par_iter().map(|g| geodetic_to_ecef(ell, g)).collect()
    }

    /// Parallel [`super::ecef_to_geodetic_batch`].
    #[must_use]
    pub fn ecef_to_geodetic_batch(ell: &Ellipsoid, points: &[Ecef]) -> Vec<Geodetic> {
        points.par_iter().map(|&p| ecef_to_geodetic(ell, p)).collect()
    }

    /// Parallel [`super::contains_batch`].
    #[must_use]
    pub fn contains_batch(mp: &MultiPolygon, points: &[Vec3]) -> Vec<bool> {
        points.par_iter().map(|&p| mp.contains(p)).collect()
    }

    /// Parallel [`super::points_in_beam`].
    ///
    /// # Errors
    /// [`GeoError::InvalidArgument`] for an invalid beam.
    pub fn points_in_beam(
        frame: &SatelliteFrame,
        boresight_body: Vec3,
        shape: &BeamShape,
        targets: &[Ecef],
    ) -> Result<Vec<bool>> {
        let beam = Beam::new(boresight_body, shape)?;
        Ok(targets.par_iter().map(|&t| beam.contains(frame, t)).collect())
    }

    /// Parallel [`super::scan_angles_batch`].
    #[must_use]
    pub fn scan_angles_batch(frame: &SatelliteFrame, targets: &[Ecef]) -> Vec<Option<ScanAngles>> {
        targets.par_iter().map(|&t| frame.scan_angles(t).ok()).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ellipsoid::WGS84;
    use crate::polygon::{lonlat_to_unit, Polygon};
    use crate::satellite::{geo_satellite_position, NadirReference};

    fn pts() -> Vec<Geodetic> {
        (0..50)
            .map(|i| Geodetic::from_degrees(-60.0 + 2.4 * f64::from(i), -170.0 + 7.0 * f64::from(i), 0.0).unwrap())
            .collect()
    }

    #[test]
    fn conversions_round_trip() {
        let g = pts();
        let e = geodetic_to_ecef_batch(&WGS84, &g);
        let back = ecef_to_geodetic_batch(&WGS84, &e);
        for (a, b) in g.iter().zip(&back) {
            assert!((a.lat_rad - b.lat_rad).abs() < 1e-12);
        }
    }

    #[test]
    fn look_angles_buffers() {
        let f = LocalFrame::new(&Geodetic::from_degrees(45.0, 0.0, 0.0).unwrap());
        let t = vec![geo_satellite_position(0.0); 3];
        let v = look_angles_batch(&f, &t);
        let mut buf = vec![v[0]; 3];
        look_angles_into(&f, &t, &mut buf).unwrap();
        assert_eq!(buf, v);
        let mut short = vec![v[0]; 2];
        assert!(matches!(
            look_angles_into(&f, &t, &mut short),
            Err(GeoError::InvalidArgument(_))
        ));
    }

    #[test]
    fn satellite_batches() {
        let f = SatelliteFrame::geo_nadir(geo_satellite_position(0.0), NadirReference::Geocentric).unwrap();
        let targets = geodetic_to_ecef_batch(&WGS84, &pts());
        let mut with_self = targets.clone();
        with_self.push(f.position);
        let s = scan_angles_batch(&f, &with_self);
        assert!(s[0].is_some());
        assert!(s.last().unwrap().is_none());
        let shape = BeamShape::Circular { half_angle_rad: 0.15 };
        let inb = points_in_beam(&f, Vec3::Z, &shape, &targets).unwrap();
        assert_eq!(inb.len(), targets.len());
        let bad = BeamShape::Circular { half_angle_rad: 0.0 };
        assert!(matches!(
            points_in_beam(&f, Vec3::Z, &bad, &targets),
            Err(GeoError::InvalidArgument(_))
        ));
    }

    #[test]
    fn contains_many() {
        let mp: MultiPolygon = Polygon::circle(Vec3::Z, 0.5, 64).unwrap().into();
        let r = contains_batch(&mp, &[lonlat_to_unit([0.0, 89.0]), lonlat_to_unit([0.0, 0.0])]);
        assert_eq!(r, vec![true, false]);
    }

    #[cfg(feature = "parallel")]
    #[test]
    fn parallel_matches_serial() {
        let g = pts();
        let e = geodetic_to_ecef_batch(&WGS84, &g);
        assert_eq!(par::geodetic_to_ecef_batch(&WGS84, &g), e);
        assert_eq!(
            par::ecef_to_geodetic_batch(&WGS84, &e),
            ecef_to_geodetic_batch(&WGS84, &e)
        );
        let f = LocalFrame::new(&g[0]);
        assert_eq!(par::look_angles_batch(&f, &e), look_angles_batch(&f, &e));
        let mp: MultiPolygon = Polygon::circle(Vec3::Z, 0.5, 64).unwrap().into();
        let u: Vec<Vec3> = g.iter().map(Geodetic::to_unit_vector).collect();
        assert_eq!(par::contains_batch(&mp, &u), contains_batch(&mp, &u));
        let sf = SatelliteFrame::geo_nadir(geo_satellite_position(0.0), NadirReference::Geocentric).unwrap();
        assert_eq!(par::scan_angles_batch(&sf, &e), scan_angles_batch(&sf, &e));
        let bad = BeamShape::Circular { half_angle_rad: -1.0 };
        assert!(matches!(
            par::points_in_beam(&sf, Vec3::Z, &bad, &e),
            Err(GeoError::InvalidArgument(_))
        ));
        let shape = BeamShape::Circular { half_angle_rad: 0.1 };
        assert_eq!(
            par::points_in_beam(&sf, Vec3::Z, &shape, &e).unwrap(),
            points_in_beam(&sf, Vec3::Z, &shape, &e).unwrap()
        );
    }
}
