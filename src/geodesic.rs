//! Distances, bearings and paths between points on the Earth.
//!
//! Two families are provided:
//! * **Spherical** (great-circle) functions: fastest, ~0.5 % error versus
//!   the ellipsoid, and naturally antimeridian- and pole-safe because they
//!   work on unit vectors.
//! * **Ellipsoidal** functions (Vincenty): millimetre accuracy on WGS-84.
//!
//! All angles are in radians, all distances in metres.

use crate::coords::Geodetic;
use crate::ellipsoid::{Ellipsoid, EARTH_MEAN_RADIUS_M};
use crate::error::{GeoError, Result};
use crate::math::{wrap_pi, wrap_two_pi, Vec3};

/// Below this projection length a point is treated as a pole of a great
/// circle (|sin| of its distance from the pole, ≈ 6 µm on the Earth).
const POLE_OF_CIRCLE_EPS: f64 = 1e-12;

/// Central angle (radians) between two points on the sphere, using the
/// unit-vector `atan2` formulation, which is accurate for all separations
/// including antipodal ones.
#[inline]
#[must_use]
pub fn central_angle(a: &Geodetic, b: &Geodetic) -> f64 {
    a.to_unit_vector().angle_to(b.to_unit_vector())
}

/// Great-circle distance in metres on a sphere of radius `radius_m`.
#[inline]
#[must_use]
pub fn great_circle_distance(a: &Geodetic, b: &Geodetic, radius_m: f64) -> f64 {
    central_angle(a, b) * radius_m
}

/// Haversine distance in metres on the mean Earth sphere.
///
/// Equivalent to [`great_circle_distance`] with [`EARTH_MEAN_RADIUS_M`];
/// kept for familiarity.
#[inline]
#[must_use]
pub fn haversine_distance(a: &Geodetic, b: &Geodetic) -> f64 {
    let dlat = b.lat_rad - a.lat_rad;
    let dlon = b.lon_rad - a.lon_rad;
    let h = (dlat * 0.5).sin().powi(2) + a.lat_rad.cos() * b.lat_rad.cos() * (dlon * 0.5).sin().powi(2);
    2.0 * EARTH_MEAN_RADIUS_M * h.sqrt().min(1.0).asin()
}

/// Initial great-circle bearing (forward azimuth) from `a` to `b`, clockwise
/// from north in `[0, 2π)`.
#[must_use]
pub fn initial_bearing(a: &Geodetic, b: &Geodetic) -> f64 {
    let dlon = b.lon_rad - a.lon_rad;
    let (sb, cb) = b.lat_rad.sin_cos();
    let (sa, ca) = a.lat_rad.sin_cos();
    let y = dlon.sin() * cb;
    let x = ca * sb - sa * cb * dlon.cos();
    wrap_two_pi(y.atan2(x))
}

/// Point reached by travelling `distance_m` along a great circle from
/// `start` with initial `bearing_rad`, on a sphere of radius `radius_m`.
/// The altitude of `start` is kept.
#[must_use]
pub fn destination(start: &Geodetic, bearing_rad: f64, distance_m: f64, radius_m: f64) -> Geodetic {
    let d = distance_m / radius_m;
    let (sd, cd) = d.sin_cos();
    let (sp, cp) = start.lat_rad.sin_cos();
    let (sb, cb) = bearing_rad.sin_cos();
    let sp2 = (sp * cd + cp * sd * cb).clamp(-1.0, 1.0);
    let lat2 = sp2.asin();
    let lon2 = start.lon_rad + (sb * sd * cp).atan2(cd - sp * sp2);
    Geodetic {
        lat_rad: lat2,
        lon_rad: wrap_pi(lon2),
        alt_m: start.alt_m,
    }
}

/// Point at fraction `t ∈ [0, 1]` along the great-circle arc from `a` to
/// `b` (spherical linear interpolation). Antimeridian- and pole-safe.
/// Altitude is linearly interpolated.
#[must_use]
pub fn interpolate(a: &Geodetic, b: &Geodetic, t: f64) -> Geodetic {
    let va = a.to_unit_vector();
    let vb = b.to_unit_vector();
    let v = slerp(va, vb, t);
    let mut g = unit_vector_to_geodetic(v);
    g.alt_m = a.alt_m + (b.alt_m - a.alt_m) * t;
    g
}

/// Midpoint of the great-circle arc between `a` and `b`.
#[must_use]
pub fn midpoint(a: &Geodetic, b: &Geodetic) -> Geodetic {
    interpolate(a, b, 0.5)
}

/// Spherical linear interpolation between two unit vectors.
///
/// For (nearly) antipodal inputs the great circle is not unique; an
/// arbitrary one is chosen.
#[must_use]
pub fn slerp(a: Vec3, b: Vec3, t: f64) -> Vec3 {
    let omega = a.angle_to(b);
    if omega < 1e-12 {
        return (a + (b - a) * t).normalize();
    }
    let so = omega.sin();
    if so < 1e-12 {
        // Antipodal: rotate around any axis orthogonal to `a`.
        let axis = a.any_orthogonal();
        return a.rotate_about(axis, omega * t);
    }
    (a * ((1.0 - t) * omega).sin() + b * (t * omega).sin()) / so
}

/// Converts a unit vector (geocentric direction) to a spherical lat/lon,
/// altitude zero.
#[inline]
#[must_use]
pub fn unit_vector_to_geodetic(v: Vec3) -> Geodetic {
    Geodetic {
        lat_rad: v.z.atan2(v.x.hypot(v.y)),
        lon_rad: v.y.atan2(v.x),
        alt_m: 0.0,
    }
}

/// Signed cross-track distance (metres, sphere) of point `p` from the great
/// circle through `a` and `b`. Positive when `p` is to the left of the path
/// direction `a → b`.
///
/// # Errors
/// [`GeoError::InvalidArgument`] when `a` and `b` coincide or are antipodal
/// (the great circle is undefined).
pub fn cross_track_distance(a: &Geodetic, b: &Geodetic, p: &Geodetic, radius_m: f64) -> Result<f64> {
    let n = great_circle_normal(a, b)?;
    Ok(cross_track_angle(n, p.to_unit_vector()) * radius_m)
}

/// Along-track distance (metres, sphere) from `a` to the foot of the
/// perpendicular dropped from `p` onto the great circle `a → b`. Negative
/// when the foot lies behind `a`.
///
/// # Errors
/// [`GeoError::InvalidArgument`] when `a` and `b` coincide or are antipodal.
pub fn along_track_distance(a: &Geodetic, b: &Geodetic, p: &Geodetic, radius_m: f64) -> Result<f64> {
    let n = great_circle_normal(a, b)?;
    Ok(along_track_angle(n, a.to_unit_vector(), p.to_unit_vector()) * radius_m)
}

/// Minimum distance (metres, sphere) from point `p` to the great-circle
/// *segment* `a → b` (not the full circle).
///
/// # Errors
/// [`GeoError::InvalidArgument`] when `a` and `b` are antipodal.
pub fn distance_to_segment(a: &Geodetic, b: &Geodetic, p: &Geodetic, radius_m: f64) -> Result<f64> {
    let (va, vb, vp) = (a.to_unit_vector(), b.to_unit_vector(), p.to_unit_vector());
    let seg = va.angle_to(vb);
    if seg < 1e-15 {
        return Ok(va.angle_to(vp) * radius_m);
    }
    let n = great_circle_normal(a, b)?;
    let along = along_track_angle(n, va, vp);
    let angle = if along > 0.0 && along < seg {
        cross_track_angle(n, vp).abs()
    } else {
        va.angle_to(vp).min(vb.angle_to(vp))
    };
    Ok(angle * radius_m)
}

/// Signed angle of `vp` from the great circle of unit normal `n`.
fn cross_track_angle(n: Vec3, vp: Vec3) -> f64 {
    n.dot(vp).clamp(-1.0, 1.0).asin()
}

/// Signed angle along the great circle of normal `n` from `va` to the foot
/// of `vp`; 0 when `vp` is a pole of the circle.
fn along_track_angle(n: Vec3, va: Vec3, vp: Vec3) -> f64 {
    let foot = vp - n * n.dot(vp);
    if foot.norm() < POLE_OF_CIRCLE_EPS {
        // Every point of the circle is equidistant from its pole.
        return 0.0;
    }
    let foot = foot.normalize();
    va.cross(foot).dot(n).atan2(va.dot(foot))
}

fn great_circle_normal(a: &Geodetic, b: &Geodetic) -> Result<Vec3> {
    let n = a.to_unit_vector().cross(b.to_unit_vector());
    // |a×b| = sin(angle): below 1e-12 the plane is numerically undefined.
    (n.norm() > 1e-12)
        .then(|| n.normalize())
        .ok_or_else(|| GeoError::InvalidArgument("great circle undefined (coincident or antipodal points)".into()))
}

/// Result of an ellipsoidal inverse geodesic problem.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeodesicInverse {
    /// Geodesic distance in metres.
    pub distance_m: f64,
    /// Forward azimuth at the first point, clockwise from north, `[0, 2π)`.
    pub azimuth1_rad: f64,
    /// Forward azimuth at the second point, clockwise from north, `[0, 2π)`.
    pub azimuth2_rad: f64,
}

/// Iteration cap for Vincenty's formulas.
const VINCENTY_MAX_ITER: usize = 200;
/// Convergence threshold on λ (inverse) or σ (direct), radians (~6 µm).
const VINCENTY_TOL: f64 = 1e-12;

/// Vincenty's series coefficients `(A, B)` for `u² = cos²α · e'²`.
fn vincenty_coefficients(ell: &Ellipsoid, cos2_alpha: f64) -> (f64, f64) {
    let u_sq = cos2_alpha * ell.ep2;
    let a = 1.0 + u_sq / 16384.0 * (4096.0 + u_sq * (-768.0 + u_sq * (320.0 - 175.0 * u_sq)));
    let b = u_sq / 1024.0 * (256.0 + u_sq * (-128.0 + u_sq * (74.0 - 47.0 * u_sq)));
    (a, b)
}

/// Angular distance on the auxiliary sphere and the terms derived from it.
#[derive(Debug, Clone, Copy)]
struct Sigma {
    sigma: f64,
    sin: f64,
    cos: f64,
    /// cos(2σₘ), σₘ being the arc midpoint's angular distance from the node.
    cos_2m: f64,
}

impl Sigma {
    /// Δσ of Vincenty's formulas.
    fn delta(&self, coeff_b: f64) -> f64 {
        let (s, c, c2m) = (self.sin, self.cos, self.cos_2m);
        coeff_b
            * s
            * (c2m
                + coeff_b / 4.0
                    * (c * (-1.0 + 2.0 * c2m * c2m)
                        - coeff_b / 6.0 * c2m * (-3.0 + 4.0 * s * s) * (-3.0 + 4.0 * c2m * c2m)))
    }

    /// Difference between longitude on the ellipsoid and on the auxiliary
    /// sphere: `(1 − C) f sinα (σ + C sinσ (cos2σₘ + C cosσ (2cos²2σₘ − 1)))`.
    fn longitude_term(&self, f: f64, sin_alpha: f64, cos2_alpha: f64) -> f64 {
        let c = f / 16.0 * cos2_alpha * (4.0 + f * (4.0 - 3.0 * cos2_alpha));
        let inner = self.cos_2m + c * self.cos * (-1.0 + 2.0 * self.cos_2m * self.cos_2m);
        (1.0 - c) * f * sin_alpha * (self.sigma + c * self.sin * inner)
    }
}

/// Vincenty's inverse formula on an ellipsoid: geodesic distance and
/// azimuths between two points (heights ignored).
///
/// # Errors
/// [`GeoError::NoConvergence`] for nearly antipodal points, where Vincenty's
/// iteration fails. Callers can fall back to [`great_circle_distance`].
pub fn vincenty_inverse(ell: &Ellipsoid, p1: &Geodetic, p2: &Geodetic) -> Result<GeodesicInverse> {
    let f = ell.f;
    let l = wrap_pi(p2.lon_rad - p1.lon_rad);
    let (su1, cu1) = ((1.0 - f) * p1.lat_rad.tan()).atan().sin_cos();
    let (su2, cu2) = ((1.0 - f) * p2.lat_rad.tan()).atan().sin_cos();

    let mut lambda = l;
    for _ in 0..VINCENTY_MAX_ITER {
        let (sl, cl) = lambda.sin_cos();
        let sin_sigma = ((cu2 * sl).powi(2) + (cu1 * su2 - su1 * cu2 * cl).powi(2)).sqrt();
        if sin_sigma == 0.0 {
            return Ok(GeodesicInverse {
                distance_m: 0.0,
                azimuth1_rad: 0.0,
                azimuth2_rad: 0.0,
            });
        }
        let cos_sigma = su1 * su2 + cu1 * cu2 * cl;
        let sin_alpha = cu1 * cu2 * sl / sin_sigma;
        let cos2_alpha = 1.0 - sin_alpha * sin_alpha;
        // On an equatorial line cos²α = 0 and cos2σₘ is taken as 0.
        let cos_2m = if cos2_alpha != 0.0 {
            cos_sigma - 2.0 * su1 * su2 / cos2_alpha
        } else {
            0.0
        };
        let sg = Sigma {
            sigma: sin_sigma.atan2(cos_sigma),
            sin: sin_sigma,
            cos: cos_sigma,
            cos_2m,
        };
        let lambda_prev = lambda;
        lambda = l + sg.longitude_term(f, sin_alpha, cos2_alpha);
        if lambda.abs() > std::f64::consts::PI + 1e-9 {
            break; // diverging: nearly antipodal
        }
        if (lambda - lambda_prev).abs() < VINCENTY_TOL {
            let (coeff_a, coeff_b) = vincenty_coefficients(ell, cos2_alpha);
            let (sl, cl) = lambda.sin_cos();
            return Ok(GeodesicInverse {
                distance_m: ell.b * coeff_a * (sg.sigma - sg.delta(coeff_b)),
                azimuth1_rad: wrap_two_pi((cu2 * sl).atan2(cu1 * su2 - su1 * cu2 * cl)),
                azimuth2_rad: wrap_two_pi((cu1 * sl).atan2(-su1 * cu2 + cu1 * su2 * cl)),
            });
        }
    }
    Err(GeoError::NoConvergence("vincenty inverse (nearly antipodal points)"))
}

/// Ellipsoidal geodesic distance in metres, falling back to the spherical
/// great-circle distance (mean radius) when Vincenty does not converge.
#[must_use]
pub fn geodesic_distance(ell: &Ellipsoid, p1: &Geodetic, p2: &Geodetic) -> f64 {
    vincenty_inverse(ell, p1, p2)
        .map(|r| r.distance_m)
        .unwrap_or_else(|_| great_circle_distance(p1, p2, EARTH_MEAN_RADIUS_M))
}

/// Vincenty's direct formula: destination reached from `start` after
/// `distance_m` metres along the geodesic with initial azimuth `azimuth_rad`.
/// Returns the destination and the forward azimuth there.
///
/// The σ iteration is a contraction (factor ≈ B ≈ 0.003 on the Earth), so it
/// always converges for finite input.
///
/// # Errors
/// [`GeoError::InvalidArgument`] for a non-finite azimuth or distance.
pub fn vincenty_direct(
    ell: &Ellipsoid,
    start: &Geodetic,
    azimuth_rad: f64,
    distance_m: f64,
) -> Result<(Geodetic, f64)> {
    if !(azimuth_rad.is_finite() && distance_m.is_finite()) {
        return Err(GeoError::InvalidArgument(format!(
            "azimuth {azimuth_rad} and distance {distance_m} must be finite"
        )));
    }
    let f = ell.f;
    let (sa1, ca1) = azimuth_rad.sin_cos();
    let tu1 = (1.0 - f) * start.lat_rad.tan();
    let cu1 = 1.0 / (1.0 + tu1 * tu1).sqrt();
    let su1 = tu1 * cu1;
    let sigma1 = tu1.atan2(ca1);
    let sin_alpha = cu1 * sa1;
    let cos2_alpha = 1.0 - sin_alpha * sin_alpha;
    let (coeff_a, coeff_b) = vincenty_coefficients(ell, cos2_alpha);

    let sigma_at = |sigma: f64| {
        let (sin, cos) = sigma.sin_cos();
        Sigma {
            sigma,
            sin,
            cos,
            cos_2m: (2.0 * sigma1 + sigma).cos(),
        }
    };
    let sigma0 = distance_m / (ell.b * coeff_a);
    let mut sigma = sigma0;
    for _ in 0..VINCENTY_MAX_ITER {
        let next = sigma0 + sigma_at(sigma).delta(coeff_b);
        let done = (next - sigma).abs() < VINCENTY_TOL;
        sigma = next;
        if done {
            break;
        }
    }
    let sg = sigma_at(sigma);
    let tmp = su1 * sg.sin - cu1 * sg.cos * ca1;
    let lat2 = (su1 * sg.cos + cu1 * sg.sin * ca1).atan2((1.0 - f) * (sin_alpha * sin_alpha + tmp * tmp).sqrt());
    let lambda = (sg.sin * sa1).atan2(cu1 * sg.cos - su1 * sg.sin * ca1);
    let l = lambda - sg.longitude_term(f, sin_alpha, cos2_alpha);
    Ok((
        Geodetic {
            lat_rad: lat2,
            lon_rad: wrap_pi(start.lon_rad + l),
            alt_m: start.alt_m,
        },
        wrap_two_pi(sin_alpha.atan2(-tmp)),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ellipsoid::WGS84;
    use approx::assert_abs_diff_eq;

    fn g(lat: f64, lon: f64) -> Geodetic {
        Geodetic::from_degrees(lat, lon, 0.0).unwrap()
    }

    #[test]
    fn haversine_matches_central_angle() {
        let a = g(48.8566, 2.3522); // Paris
        let b = g(40.7128, -74.0060); // New York
        let h = haversine_distance(&a, &b);
        let gc = great_circle_distance(&a, &b, EARTH_MEAN_RADIUS_M);
        assert_abs_diff_eq!(h, gc, epsilon = 1e-3);
        assert_abs_diff_eq!(h / 1000.0, 5837.0, epsilon = 5.0);
    }

    #[test]
    fn distance_across_antimeridian() {
        let a = g(0.0, 179.5);
        let b = g(0.0, -179.5);
        let d = great_circle_distance(&a, &b, EARTH_MEAN_RADIUS_M);
        assert_abs_diff_eq!(d, 1.0_f64.to_radians() * EARTH_MEAN_RADIUS_M, epsilon = 1e-6);
        let v = vincenty_inverse(&WGS84, &a, &b).unwrap();
        assert_abs_diff_eq!(v.distance_m, 111_319.49, epsilon = 0.01);
        assert_abs_diff_eq!(v.azimuth1_rad.to_degrees(), 90.0, epsilon = 1e-9);
    }

    #[test]
    fn bearing_and_destination_round_trip() {
        let a = g(10.0, 170.0);
        let brg = 45.0_f64.to_radians();
        let b = destination(&a, brg, 2_000_000.0, EARTH_MEAN_RADIUS_M);
        assert!(b.lon_rad < 0.0, "crossed the antimeridian");
        assert_abs_diff_eq!(
            great_circle_distance(&a, &b, EARTH_MEAN_RADIUS_M),
            2_000_000.0,
            epsilon = 1e-4
        );
        assert_abs_diff_eq!(initial_bearing(&a, &b), brg, epsilon = 1e-9);
        assert_abs_diff_eq!(initial_bearing(&g(0.0, 0.0), &g(10.0, 0.0)), 0.0, epsilon = 1e-12);
        assert_abs_diff_eq!(
            initial_bearing(&g(0.0, 0.0), &g(0.0, -10.0)).to_degrees(),
            270.0,
            epsilon = 1e-9
        );
    }

    #[test]
    fn interpolation_over_pole_and_antimeridian() {
        let a = g(80.0, 0.0);
        let b = g(80.0, 180.0);
        let m = midpoint(&a, &b);
        assert_abs_diff_eq!(m.lat_deg(), 90.0, epsilon = 1e-9);
        let a = g(0.0, 170.0);
        let b = g(0.0, -170.0);
        let m = midpoint(&a, &b);
        assert_abs_diff_eq!(m.lon_deg().abs(), 180.0, epsilon = 1e-9);
        let e = interpolate(&a, &b, 1.0);
        assert_abs_diff_eq!(e.lon_deg(), -170.0, epsilon = 1e-9);
    }

    #[test]
    fn slerp_degenerate_cases() {
        let v = slerp(Vec3::X, Vec3::X, 0.5);
        assert_abs_diff_eq!(v.x, 1.0, epsilon = 1e-15);
        let v = slerp(Vec3::X, -Vec3::X, 0.5);
        assert_abs_diff_eq!(v.norm(), 1.0, epsilon = 1e-12);
        assert_abs_diff_eq!(v.x, 0.0, epsilon = 1e-9);
    }

    #[test]
    fn cross_and_along_track() {
        let a = g(0.0, 0.0);
        let b = g(0.0, 10.0);
        let p = g(1.0, 5.0);
        let xt = cross_track_distance(&a, &b, &p, EARTH_MEAN_RADIUS_M).unwrap();
        assert_abs_diff_eq!(xt, 1.0_f64.to_radians() * EARTH_MEAN_RADIUS_M, epsilon = 1e-3);
        let at = along_track_distance(&a, &b, &p, EARTH_MEAN_RADIUS_M).unwrap();
        assert_abs_diff_eq!(at, 5.0_f64.to_radians() * EARTH_MEAN_RADIUS_M, epsilon = 50.0);
        let behind = along_track_distance(&a, &b, &g(0.0, -3.0), 1.0).unwrap();
        assert_abs_diff_eq!(behind, -3.0_f64.to_radians(), epsilon = 1e-12);
        assert_eq!(along_track_distance(&a, &b, &g(90.0, 0.0), 1.0).unwrap(), 0.0);
        assert!(matches!(
            cross_track_distance(&a, &a, &p, 1.0),
            Err(GeoError::InvalidArgument(_))
        ));
        assert!(matches!(
            along_track_distance(&a, &g(0.0, 180.0), &p, 1.0),
            Err(GeoError::InvalidArgument(_))
        ));
    }

    #[test]
    fn segment_distance() {
        let a = g(0.0, 0.0);
        let b = g(0.0, 10.0);
        let r = 1.0;
        assert_abs_diff_eq!(
            distance_to_segment(&a, &b, &g(2.0, 5.0), r).unwrap(),
            2.0_f64.to_radians(),
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(
            distance_to_segment(&a, &b, &g(0.0, 13.0), r).unwrap(),
            3.0_f64.to_radians(),
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(
            distance_to_segment(&a, &a, &g(0.0, 1.0), r).unwrap(),
            1.0_f64.to_radians(),
            epsilon = 1e-12
        );
    }

    #[test]
    fn vincenty_reference_values() {
        // Flinders Peak → Buninyong (Vincenty 1975 reference).
        let p1 = Geodetic::from_degrees(
            -(37.0 + 57.0 / 60.0 + 3.72030 / 3600.0),
            144.0 + 25.0 / 60.0 + 29.52440 / 3600.0,
            0.0,
        )
        .unwrap();
        let p2 = Geodetic::from_degrees(
            -(37.0 + 39.0 / 60.0 + 10.15610 / 3600.0),
            143.0 + 55.0 / 60.0 + 35.38390 / 3600.0,
            0.0,
        )
        .unwrap();
        let r = vincenty_inverse(&WGS84, &p1, &p2).unwrap();
        assert_abs_diff_eq!(r.distance_m, 54_972.271, epsilon = 1e-3);
        assert_abs_diff_eq!(
            r.azimuth1_rad.to_degrees(),
            306.0 + 52.0 / 60.0 + 5.37 / 3600.0,
            epsilon = 1e-4
        );
        let (dest, _) = vincenty_direct(&WGS84, &p1, r.azimuth1_rad, r.distance_m).unwrap();
        assert_abs_diff_eq!(dest.lat_rad, p2.lat_rad, epsilon = 1e-11);
        assert_abs_diff_eq!(dest.lon_rad, p2.lon_rad, epsilon = 1e-11);
    }

    #[test]
    fn vincenty_coincident_and_antipodal() {
        let a = g(12.0, 34.0);
        assert_eq!(vincenty_inverse(&WGS84, &a, &a).unwrap().distance_m, 0.0);
        let p = g(0.0, 0.0);
        let q = g(0.5, 179.7);
        assert!(matches!(
            vincenty_inverse(&WGS84, &p, &q),
            Err(GeoError::NoConvergence(_))
        ));
        let d = geodesic_distance(&WGS84, &p, &q);
        assert!(d > 19_900_000.0 && d < 20_100_000.0);
    }

    #[test]
    fn along_track_of_circle_pole_is_zero() {
        let a = g(0.0, 0.0);
        let b = g(0.0, 10.0);
        assert_eq!(
            along_track_distance(&a, &b, &g(90.0, 0.0), EARTH_MEAN_RADIUS_M).unwrap(),
            0.0
        );
        assert_eq!(
            along_track_distance(&a, &b, &g(-90.0, 0.0), EARTH_MEAN_RADIUS_M).unwrap(),
            0.0
        );
    }

    #[test]
    fn vincenty_direct_rejects_non_finite_input() {
        let a = g(10.0, 20.0);
        for (az, d) in [(f64::NAN, 1.0), (0.0, f64::INFINITY), (0.0, f64::NAN)] {
            assert!(matches!(
                vincenty_direct(&WGS84, &a, az, d),
                Err(GeoError::InvalidArgument(_))
            ));
        }
    }

    #[test]
    fn vincenty_direct_zero_distance_returns_start() {
        let a = g(-33.0, 151.0);
        let (p, az2) = vincenty_direct(&WGS84, &a, 1.0, 0.0).unwrap();
        assert_abs_diff_eq!(p.lat_rad, a.lat_rad, epsilon = 1e-15);
        assert_abs_diff_eq!(p.lon_rad, a.lon_rad, epsilon = 1e-15);
        assert_abs_diff_eq!(az2, 1.0, epsilon = 1e-12);
    }

    #[test]
    fn vincenty_equatorial_line() {
        let r = vincenty_inverse(&WGS84, &g(0.0, 0.0), &g(0.0, 90.0)).unwrap();
        // Vincenty's series truncation leaves a few micrometres here.
        assert_abs_diff_eq!(r.distance_m, WGS84.a * std::f64::consts::FRAC_PI_2, epsilon = 1e-5);
        assert_abs_diff_eq!(r.azimuth1_rad.to_degrees(), 90.0, epsilon = 1e-9);
    }

    #[test]
    fn distance_to_segment_rejects_antipodal_endpoints() {
        assert!(matches!(
            distance_to_segment(&g(0.0, 0.0), &g(0.0, 180.0), &g(1.0, 1.0), 1.0),
            Err(GeoError::InvalidArgument(_))
        ));
    }

    #[test]
    fn geodesic_distance_uses_vincenty_when_it_converges() {
        let a = g(48.8566, 2.3522);
        let b = g(40.7128, -74.0060);
        let d = geodesic_distance(&WGS84, &a, &b);
        assert_eq!(d, vincenty_inverse(&WGS84, &a, &b).unwrap().distance_m);
        assert_abs_diff_eq!(d / 1000.0, 5853.0, epsilon = 5.0);
    }
}
