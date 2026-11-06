//! Ground footprints: satellite coverage contours for an elevation mask,
//! circular / elliptical beam footprints, analytic point-in-beam tests and a
//! Gaussian main-lobe gain model.

use crate::coords::{Ecef, ecef_to_geodetic};
use crate::ellipsoid::{Ellipsoid, WGS84};
use crate::error::{GeoError, Result};
use crate::math::{Vec3, bisect_last_true};
use crate::polygon::{Orientation, Polygon};
use crate::satellite::SatelliteFrame;
use std::f64::consts::{FRAC_PI_2, TAU};

/// Number of bisection steps (2⁻⁵⁰ of the initial interval ≈ 1e-15 rad).
const BISECT_STEPS: usize = 50;

/// Shape of an antenna beam, defined in the plane tangent to the boresight
/// (so an elliptical beam is an ellipse in "tan θ" space).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BeamShape {
    /// Circular beam of the given half-angle (radians).
    Circular {
        /// Half-angle (e.g. half of the -3 dB beamwidth), radians.
        half_angle_rad: f64,
    },
    /// Elliptical beam.
    Elliptical {
        /// Half-angle along the major axis, radians.
        semi_major_rad: f64,
        /// Half-angle along the minor axis, radians.
        semi_minor_rad: f64,
        /// Angle of the major axis from the body X axis projected onto the
        /// boresight tangent plane, towards body Y, radians.
        orientation_rad: f64,
    },
}

impl BeamShape {
    fn params(&self) -> Result<(f64, f64, f64)> {
        let (a, b, o) = match *self {
            Self::Circular { half_angle_rad } => (half_angle_rad, half_angle_rad, 0.0),
            Self::Elliptical {
                semi_major_rad,
                semi_minor_rad,
                orientation_rad,
            } => (semi_major_rad, semi_minor_rad, orientation_rad),
        };
        let ok = |x: f64| x > 0.0 && x < FRAC_PI_2;
        if !(ok(a) && ok(b) && o.is_finite()) {
            return Err(GeoError::InvalidArgument("beam half-angles must be in (0, π/2)".into()));
        }
        Ok((a.tan(), b.tan(), o))
    }
}

/// A beam ready for fast repeated tests: orthonormal basis around the
/// boresight plus the ellipse axes. Build it once, test millions of points.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Beam {
    b: Vec3,
    e1: Vec3,
    e2: Vec3,
    ta: f64,
    tb: f64,
}

impl Beam {
    /// Builds a beam from a body-frame boresight and a shape.
    ///
    /// # Errors
    /// [`GeoError::InvalidArgument`] for an invalid shape or zero boresight.
    pub fn new(boresight_body: Vec3, shape: &BeamShape) -> Result<Self> {
        let (ta, tb, o) = shape.params()?;
        let b = boresight_body
            .try_normalize()
            .ok_or_else(|| GeoError::InvalidArgument("zero boresight".into()))?;
        let e1 = (Vec3::X - b * b.x)
            .try_normalize()
            .filter(|_| b.x.abs() < 1.0 - 1e-12)
            .unwrap_or_else(|| b.any_orthogonal());
        let e2 = b.cross(e1);
        let (so, co) = o.sin_cos();
        let e1r = e1 * co + e2 * so;
        Ok(Self {
            b,
            e1: e1r,
            e2: b.cross(e1r),
            ta,
            tb,
        })
    }

    /// Normalised elliptical radius of a body direction (≤ 1 inside), or
    /// `None` behind the boresight plane.
    fn radius(&self, d: Vec3) -> Option<f64> {
        let w = d.dot(self.b);
        if w <= 0.0 {
            return None;
        }
        let x = d.dot(self.e1) / (w * self.ta);
        let y = d.dot(self.e2) / (w * self.tb);
        Some(x.hypot(y))
    }

    /// Normalised radius (0 at boresight, 1 on the edge) of an ECEF target
    /// seen from `frame`, or `None` behind the antenna.
    #[inline]
    #[must_use]
    pub fn normalized_radius(&self, frame: &SatelliteFrame, target: Ecef) -> Option<f64> {
        self.radius(frame.to_body(target - frame.position))
    }

    /// True when the ECEF target is inside the beam.
    #[inline]
    #[must_use]
    pub fn contains(&self, frame: &SatelliteFrame, target: Ecef) -> bool {
        self.normalized_radius(frame, target).is_some_and(|r| r <= 1.0)
    }

    /// Direction on the beam edge at parameter `phi`, scaled by `s ∈ [0, 1]`
    /// (`s = 0` is the boresight).
    fn edge_dir(&self, phi: f64, s: f64) -> Vec3 {
        let (sp, cp) = phi.sin_cos();
        self.b + (self.e1 * (self.ta * cp) + self.e2 * (self.tb * sp)) * s
    }
}

/// True when the ECEF target lies inside the beam (analytic test, no
/// polygon).
///
/// # Errors
/// [`GeoError::InvalidArgument`] for an invalid shape or zero boresight.
pub fn point_in_beam(frame: &SatelliteFrame, boresight_body: Vec3, shape: &BeamShape, target: Ecef) -> Result<bool> {
    Ok(Beam::new(boresight_body, shape)?.contains(frame, target))
}

/// Normalised beam radius of a target (0 at boresight, 1 on the beam edge),
/// or `None` behind the antenna.
///
/// # Errors
/// [`GeoError::InvalidArgument`] for an invalid shape or zero boresight.
pub fn beam_normalized_radius(
    frame: &SatelliteFrame,
    boresight_body: Vec3,
    shape: &BeamShape,
    target: Ecef,
) -> Result<Option<f64>> {
    Ok(Beam::new(boresight_body, shape)?.normalized_radius(frame, target))
}

/// Angle (radians) between the boresight and the direction to a target.
///
/// # Errors
/// [`GeoError::InvalidArgument`] for a zero boresight or a target at the
/// platform position.
pub fn off_boresight_angle(frame: &SatelliteFrame, boresight_body: Vec3, target: Ecef) -> Result<f64> {
    let b = boresight_body
        .try_normalize()
        .ok_or_else(|| GeoError::InvalidArgument("zero boresight".into()))?;
    Ok(b.angle_to(frame.direction_to(target)?))
}

/// Gaussian main-lobe gain in dBi: `G = Gmax − 12 (θ / θ3dB)²`, where
/// `theta_3db_rad` is the full −3 dB beamwidth. Valid in the main lobe only
/// (roughly θ < θ3dB); ITU side-lobe envelopes are out of scope here.
#[must_use]
pub fn gaussian_gain_db(peak_gain_dbi: f64, theta_rad: f64, theta_3db_rad: f64) -> f64 {
    let r = theta_rad / theta_3db_rad;
    peak_gain_dbi - 12.0 * r * r
}

fn ground_point(ell: &Ellipsoid, origin: Ecef, dir_ecef: Vec3) -> Option<Vec3> {
    ell.ray_intersection(origin, dir_ecef, 0.0).map(|t| {
        let p = origin + dir_ecef * t;
        ecef_to_geodetic(ell, p).to_unit_vector()
    })
}

const BROKEN_BISECTION: &str = "bisection keeps its lower bound on the Earth";

/// Ground footprint of a beam on WGS-84, as a polygon of `n` vertices.
/// Edge rays that miss the Earth are pulled back to the limb.
///
/// # Errors
/// [`GeoError::InvalidArgument`] for `n < 3` or an invalid shape;
/// [`GeoError::NoIntersection`] when the boresight misses the Earth.
pub fn beam_footprint(frame: &SatelliteFrame, boresight_body: Vec3, shape: &BeamShape, n: usize) -> Result<Polygon> {
    if n < 3 {
        return Err(GeoError::InvalidArgument(format!(
            "footprint needs at least 3 vertices, got {n}"
        )));
    }
    let beam = Beam::new(boresight_body, shape)?;
    let hit = |phi: f64, s: f64| ground_point(&WGS84, frame.position, frame.from_body(beam.edge_dir(phi, s)));
    if hit(0.0, 0.0).is_none() {
        return Err(GeoError::NoIntersection);
    }
    let ring = (0..n)
        .map(|i| {
            let phi = TAU * i as f64 / n as f64;
            hit(phi, 1.0).unwrap_or_else(|| {
                // Edge ray misses: pull it back towards the boresight to the limb.
                let s = bisect_last_true(0.0, 1.0, BISECT_STEPS, |s| hit(phi, s).is_some());
                hit(phi, s).expect(BROKEN_BISECTION)
            })
        })
        .collect();
    Polygon::from_unit_rings(vec![ring], Orientation::SmallestArea)
}

/// Region of the Earth (WGS-84) from which the satellite is seen at or above
/// `min_elevation_rad`, as a polygon of `n` vertices.
///
/// For each azimuth around nadir, the limb is found by bisection, then the
/// nadir angle giving the requested elevation.
///
/// # Errors
/// [`GeoError::InvalidArgument`] for `n < 3`, an elevation outside
/// `[0, π/2)`, or a satellite inside the Earth.
pub fn coverage_contour(sat: Ecef, min_elevation_rad: f64, n: usize) -> Result<Polygon> {
    if n < 3 {
        return Err(GeoError::InvalidArgument(format!(
            "contour needs at least 3 vertices, got {n}"
        )));
    }
    if !(0.0..FRAC_PI_2).contains(&min_elevation_rad) {
        return Err(GeoError::InvalidArgument(format!(
            "elevation {min_elevation_rad} rad must be in [0, π/2)"
        )));
    }
    if ecef_to_geodetic(&WGS84, sat).alt_m <= 0.0 {
        return Err(GeoError::InvalidArgument(
            "satellite must be above the ellipsoid".into(),
        ));
    }
    let nadir = (-sat).normalize();
    let e1 = nadir.any_orthogonal();
    let e2 = nadir.cross(e1);
    let ring = (0..n)
        .map(|i| {
            let (sp, cp) = (TAU * i as f64 / n as f64).sin_cos();
            contour_point(sat, nadir, e1 * cp + e2 * sp, min_elevation_rad.sin())
        })
        .collect();
    Polygon::from_unit_rings(vec![ring], Orientation::SmallestArea)
}

/// Point of the coverage contour in the vertical plane spanned by `nadir`
/// and the unit vector `side`.
fn contour_point(sat: Ecef, nadir: Vec3, side: Vec3, sin_min_elevation: f64) -> Vec3 {
    let dir = |eta: f64| nadir * eta.cos() + side * eta.sin();
    let limb = bisect_last_true(0.0, FRAC_PI_2, BISECT_STEPS, |eta| {
        WGS84.ray_intersection(sat, dir(eta), 0.0).is_some()
    });
    // Elevation falls from 90° at nadir to 0° at the limb.
    let sin_elevation = |eta: f64| -> f64 {
        let d = dir(eta);
        match WGS84.ray_intersection(sat, d, 0.0) {
            Some(t) => {
                let gp = ecef_to_geodetic(&WGS84, sat + d * t);
                Ellipsoid::normal(gp.lat_rad, gp.lon_rad).dot(-d)
            }
            None => -1.0,
        }
    };
    let eta = bisect_last_true(0.0, limb, BISECT_STEPS, |eta| sin_elevation(eta) >= sin_min_elevation);
    ground_point(&WGS84, sat, dir(eta)).expect(BROKEN_BISECTION)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coords::{Geodetic, LocalFrame};
    use crate::ellipsoid::EARTH_MEAN_RADIUS_M;
    use crate::polygon::lonlat_to_unit;
    use crate::satellite::{NadirReference, ScanAngles, central_angle_for_elevation, geo_satellite_position};
    use approx::assert_abs_diff_eq;

    fn geo_frame(lon_deg: f64) -> SatelliteFrame {
        SatelliteFrame::geo_nadir(geo_satellite_position(lon_deg.to_radians()), NadirReference::Geocentric).unwrap()
    }

    #[test]
    fn coverage_contour_geo() {
        let sat = geo_satellite_position(0.0);
        let c = coverage_contour(sat, 10f64.to_radians(), 180).unwrap();
        // Every vertex sees the satellite at ~10° elevation.
        for &v in c.exterior().iter().step_by(17) {
            let ll = crate::polygon::unit_to_lonlat(v);
            let f = LocalFrame::new(&Geodetic::from_degrees(ll[1], ll[0], 0.0).unwrap());
            assert_abs_diff_eq!(f.look_angles(sat).elevation_rad.to_degrees(), 10.0, epsilon = 1e-6);
        }
        // Area close to the spherical-cap formula.
        let lambda = central_angle_for_elevation(sat.norm(), EARTH_MEAN_RADIUS_M, 10f64.to_radians());
        let cap = TAU * (1.0 - lambda.cos());
        assert_abs_diff_eq!(c.area_sr() / cap, 1.0, epsilon = 0.01);
        assert!(c.contains_lonlat([0.0, 0.0]));
        assert!(!c.contains_lonlat([100.0, 0.0]));
    }

    #[test]
    fn coverage_contour_across_antimeridian_and_errors() {
        let c = coverage_contour(geo_satellite_position(180f64.to_radians()), 0.0, 90).unwrap();
        assert!(c.contains_lonlat([-150.0, 30.0]));
        assert!(c.contains_lonlat([150.0, -30.0]));
        assert!(!c.contains_lonlat([0.0, 0.0]));
        let sat = geo_satellite_position(0.0);
        assert!(matches!(
            coverage_contour(sat, 0.1, 2),
            Err(GeoError::InvalidArgument(_))
        ));
        assert!(matches!(
            coverage_contour(sat, -0.1, 10),
            Err(GeoError::InvalidArgument(_))
        ));
        assert!(matches!(
            coverage_contour(Vec3::X, 0.1, 10),
            Err(GeoError::InvalidArgument(_))
        ));
    }

    #[test]
    fn circular_beam_footprint() {
        let f = geo_frame(0.0);
        let target = Geodetic::from_degrees(45.0, 5.0, 0.0).unwrap().to_ecef();
        let bore = f.direction_to(target).unwrap();
        let shape = BeamShape::Circular {
            half_angle_rad: 0.3f64.to_radians(),
        };
        let fp = beam_footprint(&f, bore, &shape, 72).unwrap();
        assert!(fp.contains_lonlat([5.0, 45.0]));
        assert!(!fp.contains_lonlat([5.0, 52.0]));
        // Footprint vertices are exactly on the beam edge.
        for &v in fp.exterior() {
            let ll = crate::polygon::unit_to_lonlat(v);
            let p = Geodetic::from_degrees(ll[1], ll[0], 0.0).unwrap().to_ecef();
            let r = beam_normalized_radius(&f, bore, &shape, p).unwrap().unwrap();
            assert_abs_diff_eq!(r, 1.0, epsilon = 1e-6);
        }
        assert!(point_in_beam(&f, bore, &shape, target).unwrap());
        let off = Geodetic::from_degrees(52.0, 5.0, 0.0).unwrap().to_ecef();
        assert!(!point_in_beam(&f, bore, &shape, off).unwrap());
        assert_abs_diff_eq!(off_boresight_angle(&f, bore, target).unwrap(), 0.0, epsilon = 1e-9);
    }

    #[test]
    fn elliptical_beam_orientation() {
        let f = geo_frame(0.0);
        let bore = Vec3::Z; // nadir
        let shape = BeamShape::Elliptical {
            semi_major_rad: 1f64.to_radians(),
            semi_minor_rad: 0.2f64.to_radians(),
            orientation_rad: 0.0, // major axis along body X = east
        };
        let fp = beam_footprint(&f, bore, &shape, 90).unwrap();
        // Elongated east-west.
        assert!(fp.contains_lonlat([5.0, 0.0]));
        assert!(!fp.contains_lonlat([0.0, 5.0]));
        let rotated = BeamShape::Elliptical {
            semi_major_rad: 1f64.to_radians(),
            semi_minor_rad: 0.2f64.to_radians(),
            orientation_rad: FRAC_PI_2,
        };
        let fp2 = beam_footprint(&f, bore, &rotated, 90).unwrap();
        assert!(fp2.contains_lonlat([0.0, 5.0]));
        assert!(!fp2.contains_lonlat([5.0, 0.0]));
        assert_abs_diff_eq!(fp.area_sr() / fp2.area_sr(), 1.0, epsilon = 0.01);
    }

    #[test]
    fn beam_reaching_past_the_limb() {
        let f = geo_frame(0.0);
        // Boresight near the limb with a wide beam: part of the cone misses.
        let bore = ScanAngles::from_theta_phi(8.0f64.to_radians(), 0.0).direction();
        let shape = BeamShape::Circular {
            half_angle_rad: 2f64.to_radians(),
        };
        let fp = beam_footprint(&f, bore, &shape, 64).unwrap();
        assert!(fp.area_sr() > 0.0);
        assert!(fp.contains(lonlat_to_unit([60.0, 0.0])));
    }

    #[test]
    fn beam_errors() {
        let f = geo_frame(0.0);
        let ok = BeamShape::Circular { half_angle_rad: 0.01 };
        assert!(matches!(
            beam_footprint(&f, Vec3::Z, &ok, 2),
            Err(GeoError::InvalidArgument(_))
        ));
        assert!(matches!(
            beam_footprint(&f, -Vec3::Z, &ok, 8),
            Err(GeoError::NoIntersection)
        ));
        let bad = BeamShape::Circular { half_angle_rad: 2.0 };
        assert!(matches!(
            beam_footprint(&f, Vec3::Z, &bad, 8),
            Err(GeoError::InvalidArgument(_))
        ));
        assert!(matches!(
            point_in_beam(&f, Vec3::ZERO, &ok, Vec3::X),
            Err(GeoError::InvalidArgument(_))
        ));
        assert!(matches!(
            off_boresight_angle(&f, Vec3::ZERO, Vec3::X),
            Err(GeoError::InvalidArgument(_))
        ));
        // Target behind the antenna.
        let behind = f.position * 2.0;
        assert_eq!(beam_normalized_radius(&f, Vec3::Z, &ok, behind).unwrap(), None);
        // Boresight along body X uses a fallback basis.
        assert!(point_in_beam(&f, Vec3::X, &ok, f.position + f.x_axis() * 1e6).unwrap());
    }

    #[test]
    fn gaussian_gain() {
        assert_abs_diff_eq!(gaussian_gain_db(40.0, 0.0, 0.01), 40.0);
        assert_abs_diff_eq!(gaussian_gain_db(40.0, 0.005, 0.01), 37.0);
    }

    #[test]
    fn normalized_radius_rejects_invalid_beam() {
        let f = geo_frame(0.0);
        let bad = BeamShape::Elliptical {
            semi_major_rad: 0.1,
            semi_minor_rad: 0.0,
            orientation_rad: 0.0,
        };
        assert!(matches!(
            beam_normalized_radius(&f, Vec3::Z, &bad, Vec3::X),
            Err(GeoError::InvalidArgument(_))
        ));
        let nan = BeamShape::Elliptical {
            semi_major_rad: 0.1,
            semi_minor_rad: 0.1,
            orientation_rad: f64::NAN,
        };
        assert!(matches!(Beam::new(Vec3::Z, &nan), Err(GeoError::InvalidArgument(_))));
    }

    #[test]
    fn off_boresight_angle_rejects_target_at_platform() {
        let f = geo_frame(0.0);
        assert!(matches!(
            off_boresight_angle(&f, Vec3::Z, f.position),
            Err(GeoError::InvalidArgument(_))
        ));
    }

    #[test]
    fn beam_reports_radius_one_on_edge_and_zero_on_boresight() {
        let f = geo_frame(0.0);
        let shape = BeamShape::Circular { half_angle_rad: 0.01 };
        let beam = Beam::new(Vec3::Z, &shape).unwrap();
        let bore_target = f.ground_intersection(Vec3::Z, &WGS84, 0.0).unwrap();
        assert_abs_diff_eq!(beam.normalized_radius(&f, bore_target).unwrap(), 0.0, epsilon = 1e-12);
        let edge = ScanAngles::from_theta_phi(0.01, 1.0).direction();
        let edge_target = f.ground_intersection(edge, &WGS84, 0.0).unwrap();
        assert_abs_diff_eq!(beam.normalized_radius(&f, edge_target).unwrap(), 1.0, epsilon = 1e-9);
    }
}
