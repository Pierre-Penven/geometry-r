//! Low-precision time, Earth-rotation and Sun models, with sun-outage and
//! eclipse predictions.
//!
//! Accuracy is about 0.01° for the Sun direction between 1950 and 2050.
//! Precession, nutation, polar motion and UT1−UTC are ignored (the error is
//! well below the Sun's 0.27° angular radius), so these helpers are suited to
//! outage and eclipse prediction, not to precise astrometry.

use crate::coords::Ecef;
use crate::ellipsoid::EARTH_MEAN_RADIUS_M;
use crate::error::{GeoError, Result};
use crate::math::{Mat3, Vec3, bisect_last_true, bisect_steps_for, wrap_two_pi};

/// Astronomical unit in metres.
pub const AU_M: f64 = 149_597_870_700.0;
/// Mean solar radius in metres.
pub const SUN_RADIUS_M: f64 = 695_700_000.0;
/// Julian date of the Unix epoch (1970-01-01T00:00:00Z).
pub const JD_UNIX_EPOCH: f64 = 2_440_587.5;
/// Julian date of J2000.0.
pub const JD_J2000: f64 = 2_451_545.0;

/// Accuracy (seconds) to which window edges are refined.
const WINDOW_EDGE_TOL_S: f64 = 0.1;

/// Unix time (seconds, UTC) → Julian date.
#[inline]
#[must_use]
pub fn unix_to_julian_date(unix_s: f64) -> f64 {
    unix_s / 86_400.0 + JD_UNIX_EPOCH
}

/// Greenwich mean sidereal time (radians, `[0, 2π)`) for a Julian date
/// (UT1 ≈ UTC), IAU 1982 expression.
#[must_use]
pub fn gmst_rad(jd: f64) -> f64 {
    let d = jd - JD_J2000;
    let t = d / 36_525.0;
    let deg = 280.460_618_37 + 360.985_647_366_29 * d + 0.000_387_933 * t * t - t * t * t / 38_710_000.0;
    wrap_two_pi(deg.to_radians())
}

/// Rotation matrix taking inertial (true-of-date equatorial) vectors to ECEF
/// at the given Unix time.
#[must_use]
pub fn eci_to_ecef_matrix(unix_s: f64) -> Mat3 {
    Mat3::rot_z(-gmst_rad(unix_to_julian_date(unix_s)))
}

/// Inertial → ECEF.
#[must_use]
pub fn eci_to_ecef(v: Vec3, unix_s: f64) -> Vec3 {
    eci_to_ecef_matrix(unix_s).mul_vec(v)
}

/// ECEF → inertial.
#[must_use]
pub fn ecef_to_eci(v: Vec3, unix_s: f64) -> Vec3 {
    eci_to_ecef_matrix(unix_s).mul_vec_t(v)
}

/// Sun position (metres) in the inertial equatorial frame of date.
#[must_use]
pub fn sun_position_eci(jd: f64) -> Vec3 {
    let n = jd - JD_J2000;
    let l = (280.460 + 0.985_647_4 * n).to_radians();
    let g = (357.528 + 0.985_600_3 * n).to_radians();
    let lambda = l + (1.915f64.to_radians()) * g.sin() + (0.020f64.to_radians()) * (2.0 * g).sin();
    let eps = (23.439 - 0.000_000_4 * n).to_radians();
    let r = (1.000_14 - 0.016_71 * g.cos() - 0.000_14 * (2.0 * g).cos()) * AU_M;
    let (sl, cl) = lambda.sin_cos();
    let (se, ce) = eps.sin_cos();
    Vec3::new(r * cl, r * ce * sl, r * se * sl)
}

/// Sun position in ECEF (metres) at a Unix time.
#[must_use]
pub fn sun_position_ecef(unix_s: f64) -> Ecef {
    eci_to_ecef(sun_position_eci(unix_to_julian_date(unix_s)), unix_s)
}

/// Angle (radians) between the satellite and the Sun seen from a station.
#[must_use]
pub fn sun_satellite_separation(station: Ecef, sat: Ecef, sun: Ecef) -> f64 {
    (sat - station).angle_to(sun - station)
}

/// True when the Sun is within `threshold_rad` of the satellite direction
/// seen from the station (sun outage / sun transit).
#[must_use]
pub fn is_sun_outage(station: Ecef, sat: Ecef, unix_s: f64, threshold_rad: f64) -> bool {
    sun_satellite_separation(station, sat, sun_position_ecef(unix_s)) <= threshold_rad
}

/// Predicts sun-outage windows `[start, end]` (Unix seconds) between
/// `start_s` and `end_s`.
///
/// `sat_at` returns the satellite ECEF position at a time (constant for a GEO
/// satellite). The interval is scanned every `step_s` seconds and each edge
/// is refined by bisection to 0.1 s. `step_s` must be shorter than the
/// shortest outage to be detected (a few minutes is typical).
///
/// # Errors
/// [`GeoError::InvalidArgument`] for a non-positive step or an empty span.
pub fn sun_outage_windows<F>(
    station: Ecef,
    sat_at: F,
    start_s: f64,
    end_s: f64,
    step_s: f64,
    threshold_rad: f64,
) -> Result<Vec<(f64, f64)>>
where
    F: Fn(f64) -> Ecef,
{
    find_windows(start_s, end_s, step_s, |t| {
        sun_satellite_separation(station, sat_at(t), sun_position_ecef(t)) <= threshold_rad
    })
}

/// Generic window finder on a boolean predicate of time.
pub(crate) fn find_windows<P: Fn(f64) -> bool>(
    start_s: f64,
    end_s: f64,
    step_s: f64,
    pred: P,
) -> Result<Vec<(f64, f64)>> {
    if step_s.is_nan() || step_s <= 0.0 || end_s.is_nan() || start_s.is_nan() || end_s <= start_s {
        return Err(GeoError::InvalidArgument("need step > 0 and end > start".into()));
    }
    let refine = |lo: f64, hi: f64, lo_val: bool| {
        let steps = bisect_steps_for(hi - lo, WINDOW_EDGE_TOL_S);
        bisect_last_true(lo, hi, steps, |t| pred(t) == lo_val)
    };
    let mut out = Vec::new();
    let mut t = start_s;
    let mut prev = pred(t);
    let mut open = prev.then_some(start_s);
    while t < end_s {
        let next_t = (t + step_s).min(end_s);
        let cur = pred(next_t);
        if cur != prev {
            let edge = refine(t, next_t, prev);
            if cur {
                open = Some(edge);
            } else {
                // A falling edge always follows the start or a rising edge.
                let start = open.take().expect("falling edge without an open window");
                out.push((start, edge));
            }
        }
        prev = cur;
        t = next_t;
    }
    if let Some(s) = open {
        out.push((s, end_s));
    }
    Ok(out)
}

/// Illumination state of a satellite.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EclipseState {
    /// Full sunlight.
    Sunlit,
    /// Partially shadowed by the Earth.
    Penumbra,
    /// Fully shadowed by the Earth.
    Umbra,
}

/// Conical-shadow eclipse model (spherical Earth of mean radius, finite Sun).
#[must_use]
pub fn eclipse_state(sat: Ecef, sun: Ecef) -> EclipseState {
    let to_sun = sun - sat;
    let r_sat = sat.norm();
    if r_sat <= EARTH_MEAN_RADIUS_M {
        return EclipseState::Umbra;
    }
    let a_sun = (SUN_RADIUS_M / to_sun.norm()).min(1.0).asin();
    let a_earth = (EARTH_MEAN_RADIUS_M / r_sat).min(1.0).asin();
    let sep = to_sun.angle_to(-sat);
    if sep >= a_sun + a_earth {
        EclipseState::Sunlit
    } else if sep <= a_earth - a_sun {
        EclipseState::Umbra
    } else {
        EclipseState::Penumbra
    }
}

/// Cylindrical-shadow test (faster, slightly pessimistic at the shadow
/// edge): true when the satellite is behind the Earth relative to the Sun
/// and within one Earth radius of the Earth–Sun axis.
#[must_use]
pub fn in_cylindrical_shadow(sat: Ecef, sun: Ecef) -> bool {
    let s = sun.normalize();
    let along = sat.dot(s);
    along < 0.0 && (sat - s * along).norm() < EARTH_MEAN_RADIUS_M
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coords::Geodetic;
    use crate::satellite::geo_satellite_position;
    use approx::assert_abs_diff_eq;

    // 2024-03-20T00:00:00Z and 2024-06-20T12:00:00Z.
    const MAR20: f64 = 1_710_892_800.0;
    const JUN20_NOON: f64 = 1_718_884_800.0;

    #[test]
    fn time_conversions() {
        assert_eq!(unix_to_julian_date(0.0), JD_UNIX_EPOCH);
        assert_abs_diff_eq!(gmst_rad(JD_J2000).to_degrees(), 280.460_618_37, epsilon = 1e-9);
        let v = Vec3::new(1.0, 2.0, 3.0);
        let back = ecef_to_eci(eci_to_ecef(v, MAR20), MAR20);
        assert_abs_diff_eq!(back.distance(v), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn sun_declination_and_longitude() {
        // Equinox (2024-03-20 03:06 UTC): declination ≈ 0.
        let s = sun_position_ecef(MAR20 + 3.1 * 3600.0);
        let dec = s.z.atan2(s.x.hypot(s.y)).to_degrees();
        assert_abs_diff_eq!(dec, 0.0, epsilon = 0.05);
        // June solstice: declination ≈ 23.44°, sub-solar longitude ≈ 0 at noon UTC.
        let s = sun_position_ecef(JUN20_NOON);
        let dec = s.z.atan2(s.x.hypot(s.y)).to_degrees();
        assert_abs_diff_eq!(dec, 23.44, epsilon = 0.05);
        let lon = s.y.atan2(s.x).to_degrees();
        assert_abs_diff_eq!(lon, 0.4, epsilon = 0.5); // equation of time ≈ −1.7 min
        assert_abs_diff_eq!(s.norm() / AU_M, 1.016, epsilon = 0.002);
    }

    #[test]
    fn sun_outage_at_equinox() {
        let station = Geodetic::from_degrees(0.0, 0.0, 0.0).unwrap().to_ecef();
        let sat = geo_satellite_position(0.0);
        let w = sun_outage_windows(station, |_| sat, MAR20, MAR20 + 86_400.0, 60.0, 1f64.to_radians()).unwrap();
        assert_eq!(w.len(), 1, "{w:?}");
        let mid = 0.5 * (w[0].0 + w[0].1);
        // Around local noon (12:07 UTC with the equation of time).
        assert_abs_diff_eq!((mid - MAR20) / 3600.0, 12.12, epsilon = 0.1);
        assert!(is_sun_outage(station, sat, mid, 1f64.to_radians()));
        assert!(!is_sun_outage(station, sat, MAR20, 1f64.to_radians()));
        assert!(matches!(
            sun_outage_windows(station, |_| sat, MAR20, MAR20, 60.0, 0.01),
            Err(GeoError::InvalidArgument(_))
        ));
        assert!(matches!(
            sun_outage_windows(station, |_| sat, MAR20, MAR20 + 10.0, 0.0, 0.01),
            Err(GeoError::InvalidArgument(_))
        ));
    }

    #[test]
    fn window_finder_edges() {
        // Predicate true at start and end: windows are clipped to the span.
        let w = find_windows(0.0, 100.0, 7.0, |t| !(40.0..60.0).contains(&t)).unwrap();
        assert_eq!(w.len(), 2);
        assert_eq!(w[0].0, 0.0);
        assert_abs_diff_eq!(w[0].1, 40.0, epsilon = 0.1);
        assert_abs_diff_eq!(w[1].0, 60.0, epsilon = 0.1);
        assert_eq!(w[1].1, 100.0);
    }

    #[test]
    fn geo_eclipse_at_equinox() {
        let sat = geo_satellite_position(0.0);
        // Local midnight at the satellite's longitude: in umbra.
        assert_eq!(eclipse_state(sat, sun_position_ecef(MAR20)), EclipseState::Umbra);
        assert!(in_cylindrical_shadow(sat, sun_position_ecef(MAR20)));
        // Noon: sunlit.
        assert_eq!(
            eclipse_state(sat, sun_position_ecef(MAR20 + 43_200.0)),
            EclipseState::Sunlit
        );
        assert!(!in_cylindrical_shadow(sat, sun_position_ecef(MAR20 + 43_200.0)));
        // Solstice midnight: GEO is never eclipsed.
        assert_eq!(
            eclipse_state(sat, sun_position_ecef(JUN20_NOON - 43_200.0)),
            EclipseState::Sunlit
        );
        // Inside the Earth: umbra by definition.
        assert_eq!(eclipse_state(Vec3::X, sun_position_ecef(MAR20)), EclipseState::Umbra);
        // Penumbra exists at the shadow edge: scan around midnight.
        let mut saw_penumbra = false;
        for i in 0..400 {
            let t = MAR20 - 6000.0 + f64::from(i) * 30.0;
            if eclipse_state(sat, sun_position_ecef(t)) == EclipseState::Penumbra {
                saw_penumbra = true;
            }
        }
        assert!(saw_penumbra);
    }
}
