//! Ground-station-centred geometry: look angles, visibility (with horizon
//! masks), line of sight, topocentric separations, separation from the
//! geostationary arc, and polarization skew.

use crate::coords::{Ecef, Geodetic, LocalFrame, LookAngles};
use crate::ellipsoid::WGS84;
use crate::error::{GeoError, Result};
use crate::math::{Vec3, wrap_pi, wrap_two_pi};
use crate::satellite::geo_satellite_position;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

/// Look angles from a geodetic station to an ECEF target.
///
/// For many targets seen from the same station, build a [`LocalFrame`] once
/// and call [`LocalFrame::look_angles`].
#[must_use]
pub fn look_angles(station: &Geodetic, target: Ecef) -> LookAngles {
    LocalFrame::new(station).look_angles(target)
}

/// True when `target` is at or above `min_elevation_rad` seen from the frame
/// origin.
#[inline]
#[must_use]
pub fn is_visible(frame: &LocalFrame, target: Ecef, min_elevation_rad: f64) -> bool {
    frame.look_angles(target).elevation_rad >= min_elevation_rad
}

/// A per-azimuth minimum-elevation mask (terrain, buildings, regulatory
/// mask). Values are linearly interpolated between samples, wrapping at 360°.
#[derive(Debug, Clone, PartialEq)]
pub struct HorizonMask {
    /// `(azimuth_rad, min_elevation_rad)` sorted by azimuth in `[0, 2π)`.
    points: Vec<(f64, f64)>,
}

impl HorizonMask {
    /// Builds a mask from `(azimuth_rad, min_elevation_rad)` samples, in any
    /// order. Azimuths are wrapped to `[0, 2π)`.
    ///
    /// # Errors
    /// [`GeoError::InvalidArgument`] if no sample is given, a value is not
    /// finite, or two samples share an azimuth.
    pub fn new(samples: &[(f64, f64)]) -> Result<Self> {
        if samples.is_empty() {
            return Err(GeoError::InvalidArgument("empty horizon mask".into()));
        }
        let mut points = Vec::with_capacity(samples.len());
        for &(az, el) in samples {
            if !(az.is_finite() && el.is_finite()) {
                return Err(GeoError::InvalidArgument("non-finite horizon mask sample".into()));
            }
            points.push((wrap_two_pi(az), el));
        }
        points.sort_by(|a, b| a.0.total_cmp(&b.0));
        if points.windows(2).any(|w| w[0].0 == w[1].0) {
            return Err(GeoError::InvalidArgument("duplicate azimuth in horizon mask".into()));
        }
        Ok(Self { points })
    }

    /// Minimum elevation (radians) at an azimuth (radians).
    #[must_use]
    pub fn min_elevation(&self, azimuth_rad: f64) -> f64 {
        let az = wrap_two_pi(azimuth_rad);
        let n = self.points.len();
        if n == 1 {
            return self.points[0].1;
        }
        let idx = self.points.partition_point(|p| p.0 <= az);
        let (a, b) = if idx == 0 || idx == n {
            // Between the last and the first sample, across 0/360°.
            let a = self.points[n - 1];
            let b = self.points[0];
            (a, (b.0 + TAU, b.1))
        } else {
            (self.points[idx - 1], self.points[idx])
        };
        let az = if az < a.0 { az + TAU } else { az };
        let t = (az - a.0) / (b.0 - a.0);
        a.1 + (b.1 - a.1) * t
    }

    /// True when `look` clears the mask.
    #[must_use]
    pub fn is_clear(&self, look: &LookAngles) -> bool {
        look.elevation_rad >= self.min_elevation(look.azimuth_rad)
    }
}

/// True when the straight line between two ECEF points clears the WGS-84
/// ellipsoid raised by `grazing_height_m` (e.g. 80 km to stay above the
/// dense atmosphere on inter-satellite links).
#[must_use]
pub fn line_of_sight(a: Ecef, b: Ecef, grazing_height_m: f64) -> bool {
    !WGS84.segment_blocked(a, b, grazing_height_m)
}

/// Angle (radians) between two targets as seen from an observer.
#[must_use]
pub fn topocentric_separation(observer: Ecef, a: Ecef, b: Ecef) -> f64 {
    (a - observer).angle_to(b - observer)
}

/// Result of a GSO-arc separation search.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GsoSeparation {
    /// Smallest angle (radians) between the pointing direction and the
    /// visible geostationary arc.
    pub angle_rad: f64,
    /// Longitude (radians) of the closest geostationary position.
    pub geo_lon_rad: f64,
}

/// Minimum angle between the pointing direction `station → target` and the
/// part of the geostationary arc visible above `min_elevation_rad`.
///
/// Returns `None` when no part of the arc is visible (polar stations).
/// The search samples the arc every 0.25° and refines with golden-section
/// search, giving ~1e-9 rad accuracy.
#[must_use]
pub fn gso_arc_separation(station: &LocalFrame, target: Ecef, min_elevation_rad: f64) -> Option<GsoSeparation> {
    let dir = target - station.origin;
    // Separation to the GEO slot at `lon`, infinite when the slot is masked.
    let separation = |lon: f64| -> f64 {
        let g = geo_satellite_position(lon);
        if station.look_angles(g).elevation_rad >= min_elevation_rad {
            dir.angle_to(g - station.origin)
        } else {
            f64::INFINITY
        }
    };
    let (coarse_sep, coarse_lon) = coarse_arc_scan(station.geodetic.lon_rad, &separation)?;
    let fine_lon = golden_section_min(
        coarse_lon - GSO_SCAN_STEP_RAD,
        coarse_lon + GSO_SCAN_STEP_RAD,
        &separation,
    );
    let fine_sep = separation(fine_lon);
    let (angle_rad, lon) = if fine_sep < coarse_sep {
        (fine_sep, fine_lon)
    } else {
        (coarse_sep, coarse_lon)
    };
    Some(GsoSeparation {
        angle_rad,
        geo_lon_rad: wrap_pi(lon),
    })
}

/// Sampling step of the coarse GSO-arc scan.
const GSO_SCAN_STEP_RAD: f64 = 0.25 * PI / 180.0;
/// Golden-section iterations (0.618⁶⁰ of the bracket ≈ 1e-13).
const GOLDEN_ITERATIONS: usize = 60;

/// Best `(separation, longitude)` on a grid spanning ±90° of longitude
/// around the station, or `None` when every slot is masked.
fn coarse_arc_scan(station_lon: f64, separation: &impl Fn(f64) -> f64) -> Option<(f64, f64)> {
    let n = (FRAC_PI_2 / GSO_SCAN_STEP_RAD).round() as i32;
    (-n..=n)
        .map(|i| station_lon + f64::from(i) * GSO_SCAN_STEP_RAD)
        .map(|lon| (separation(lon), lon))
        .filter(|(s, _)| s.is_finite())
        .min_by(|a, b| a.0.total_cmp(&b.0))
}

/// Golden-section search for the minimum of a unimodal `f` on `[lo, hi]`.
fn golden_section_min(mut lo: f64, mut hi: f64, f: &impl Fn(f64) -> f64) -> f64 {
    let gr = (5f64.sqrt() - 1.0) / 2.0;
    let mut c = hi - gr * (hi - lo);
    let mut d = lo + gr * (hi - lo);
    let (mut fc, mut fd) = (f(c), f(d));
    for _ in 0..GOLDEN_ITERATIONS {
        if fc < fd {
            (hi, d, fd) = (d, c, fc);
            c = hi - gr * (hi - lo);
            fc = f(c);
        } else {
            (lo, c, fc) = (c, d, fd);
            d = lo + gr * (hi - lo);
            fd = f(d);
        }
    }
    0.5 * (lo + hi)
}

/// Polarization skew angle (radians) for a linearly polarised link.
///
/// `sat_pol_ref` is the satellite's polarization reference direction in
/// ECEF (for a GEO satellite whose "vertical" polarization is aligned with
/// the Earth's axis, pass `Vec3::Z`). The result is the rotation from the
/// station's local vertical (projected onto the wavefront plane) to the
/// satellite polarization vector, **positive clockwise as seen by an
/// observer at the station facing the satellite**, wrapped to
/// `(-π/2, π/2]` (linear polarization is 180°-periodic).
///
/// # Errors
/// [`GeoError::InvalidArgument`] if the reference vector or the local
/// vertical is parallel to the line of sight (e.g. satellite at zenith with
/// an undefined reference).
pub fn polarization_skew(station: &LocalFrame, sat: Ecef, sat_pol_ref: Vec3) -> Result<f64> {
    let k = (sat - station.origin)
        .try_normalize()
        .ok_or_else(|| GeoError::InvalidArgument("satellite coincides with station".into()))?;
    let project = |v: Vec3| (v - k * v.dot(k)).try_normalize().filter(|p| p.is_finite());
    let p = project(sat_pol_ref)
        .ok_or_else(|| GeoError::InvalidArgument("polarization reference along line of sight".into()))?;
    let u = project(station.up())
        .ok_or_else(|| GeoError::InvalidArgument("satellite at zenith: vertical undefined".into()))?;
    let mut a = k.dot(u.cross(p)).atan2(u.dot(p));
    if a > FRAC_PI_2 {
        a -= PI;
    } else if a <= -FRAC_PI_2 {
        a += PI;
    }
    Ok(a)
}

/// Classic closed-form GEO polarization skew `atan(sin Δλ / tan φ)` with
/// `Δλ = λ_sat − λ_station`, in `(-π/2, π/2)`. Sign conventions differ
/// between operators; [`polarization_skew`] is the reference.
#[must_use]
pub fn geo_polarization_skew(station: &Geodetic, sat_lon_rad: f64) -> f64 {
    let dl = wrap_pi(sat_lon_rad - station.lon_rad);
    (dl.sin() / station.lat_rad.tan()).atan()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    fn st(lat: f64, lon: f64) -> Geodetic {
        Geodetic::from_degrees(lat, lon, 0.0).unwrap()
    }

    #[test]
    fn geo_look_angles_reference() {
        // Station on the equator under the satellite: zenith.
        let la = look_angles(&st(0.0, 13.0), geo_satellite_position(13f64.to_radians()));
        assert_abs_diff_eq!(la.elevation_rad.to_degrees(), 90.0, epsilon = 1e-6);
        assert_abs_diff_eq!(la.range_m, 35_786_035.0, epsilon = 1.0);
        // Paris → 13°E: elevation ~33.0°, azimuth ~166° (south-southeast).
        let la = look_angles(&st(48.85, 2.35), geo_satellite_position(13f64.to_radians()));
        assert_abs_diff_eq!(la.elevation_rad.to_degrees(), 33.0, epsilon = 0.1);
        assert_abs_diff_eq!(la.azimuth_rad.to_degrees(), 165.9, epsilon = 0.5);
        let f = LocalFrame::new(&st(48.85, 2.35));
        assert!(is_visible(&f, geo_satellite_position(13f64.to_radians()), 0.1));
        assert!(!is_visible(&f, geo_satellite_position(-150f64.to_radians()), 0.0));
    }

    #[test]
    fn horizon_mask_interpolation() {
        let m = HorizonMask::new(&[(0.0, 0.1), (FRAC_PI_2, 0.3), (PI, 0.1)]).unwrap();
        assert_abs_diff_eq!(m.min_elevation(FRAC_PI_2 / 2.0), 0.2, epsilon = 1e-12);
        // Wrap-around between 180° and 360°/0°.
        assert_abs_diff_eq!(m.min_elevation(1.5 * PI), 0.1, epsilon = 1e-12);
        assert_abs_diff_eq!(m.min_elevation(-0.0), 0.1, epsilon = 1e-12);
        let single = HorizonMask::new(&[(1.0, 0.05)]).unwrap();
        assert_eq!(single.min_elevation(3.0), 0.05);
        let la = LookAngles {
            azimuth_rad: FRAC_PI_2,
            elevation_rad: 0.25,
            range_m: 1.0,
        };
        assert!(!m.is_clear(&la));
        assert!(matches!(HorizonMask::new(&[]), Err(GeoError::InvalidArgument(_))));
        assert!(matches!(
            HorizonMask::new(&[(0.0, 0.0), (TAU, 0.1)]),
            Err(GeoError::InvalidArgument(_))
        ));
        assert!(matches!(
            HorizonMask::new(&[(f64::NAN, 0.0)]),
            Err(GeoError::InvalidArgument(_))
        ));
    }

    #[test]
    fn line_of_sight_cases() {
        let a = Geodetic::from_degrees(0.0, 0.0, 550e3).unwrap().to_ecef();
        let b = Geodetic::from_degrees(0.0, 20.0, 550e3).unwrap().to_ecef();
        let c = Geodetic::from_degrees(0.0, 90.0, 550e3).unwrap().to_ecef();
        assert!(line_of_sight(a, b, 80e3));
        assert!(!line_of_sight(a, c, 0.0));
    }

    #[test]
    fn separation_between_geo_satellites() {
        let obs = st(0.0, 0.0).to_ecef();
        let s = topocentric_separation(
            obs,
            geo_satellite_position(0.0),
            geo_satellite_position(2f64.to_radians()),
        );
        // Topocentric spacing is larger than geocentric spacing (~2.36° for 2°).
        assert_abs_diff_eq!(s.to_degrees(), 2.0 * 42_164.0 / 35_786.0, epsilon = 0.01);
    }

    #[test]
    fn gso_arc_separation_cases() {
        let f = LocalFrame::new(&st(45.0, 10.0));
        // Pointing exactly at a GEO slot: zero separation.
        let r = gso_arc_separation(&f, geo_satellite_position(20f64.to_radians()), 0.0).unwrap();
        assert_abs_diff_eq!(r.angle_rad, 0.0, epsilon = 1e-7);
        assert_abs_diff_eq!(r.geo_lon_rad.to_degrees(), 20.0, epsilon = 1e-4);
        // Pointing at zenith: separation is ~ 90° − elevation of arc's highest point.
        let zen = f.from_enu(Vec3::new(0.0, 0.0, 1e6));
        let r = gso_arc_separation(&f, zen, 0.0).unwrap();
        let max_el = look_angles(&st(45.0, 10.0), geo_satellite_position(10f64.to_radians())).elevation_rad;
        assert_abs_diff_eq!(r.angle_rad, FRAC_PI_2 - max_el, epsilon = 1e-6);
        // Polar station sees no arc.
        let polar = LocalFrame::new(&st(89.0, 0.0));
        assert!(gso_arc_separation(&polar, zen, 0.0).is_none());
    }

    #[test]
    fn polarization_skew_matches_formula() {
        for &(lat, lon, sat) in &[
            (45.0, 0.0, 0.0),
            (45.0, 10.0, 0.0),
            (30.0, -20.0, 5.0),
            (-35.0, 20.0, 0.0),
        ] {
            let g = st(lat, lon);
            let f = LocalFrame::new(&g);
            let v = polarization_skew(&f, geo_satellite_position(f64::to_radians(sat)), Vec3::Z).unwrap();
            let c = geo_polarization_skew(&g, f64::to_radians(sat));
            // Same magnitude (the closed form is a spherical approximation).
            assert_abs_diff_eq!(v.abs().to_degrees(), c.abs().to_degrees(), epsilon = 1.0);
        }
        let g = st(45.0, 0.0);
        let f = LocalFrame::new(&g);
        let v = polarization_skew(&f, geo_satellite_position(0.0), Vec3::Z).unwrap();
        assert_abs_diff_eq!(v, 0.0, epsilon = 1e-9);
        // Mirror symmetry: east and west give opposite signs.
        let e = polarization_skew(&f, geo_satellite_position(0.2), Vec3::Z).unwrap();
        let w = polarization_skew(&f, geo_satellite_position(-0.2), Vec3::Z).unwrap();
        assert_abs_diff_eq!(e, -w, epsilon = 1e-9);
        assert!(e.abs() > 0.05);
        let k = geo_satellite_position(0.0) - f.origin;
        assert!(matches!(
            polarization_skew(&f, geo_satellite_position(0.0), k),
            Err(GeoError::InvalidArgument(_))
        ));
        let eq = LocalFrame::new(&st(0.0, 0.0));
        assert!(matches!(
            polarization_skew(&eq, geo_satellite_position(0.0), Vec3::Z),
            Err(GeoError::InvalidArgument(_))
        ));
        assert!(matches!(
            polarization_skew(&f, f.origin, Vec3::Z),
            Err(GeoError::InvalidArgument(_))
        ));
    }

    #[test]
    fn polarization_skew_is_180_periodic_in_reference_vector() {
        let f = LocalFrame::new(&st(40.0, 0.0));
        for lon in [-0.4, -0.1, 0.1, 0.4] {
            let sat = geo_satellite_position(lon);
            let up = polarization_skew(&f, sat, Vec3::Z).unwrap();
            let down = polarization_skew(&f, sat, -Vec3::Z).unwrap();
            assert_abs_diff_eq!(up, down, epsilon = 1e-12);
            assert!(up > -FRAC_PI_2 && up <= FRAC_PI_2);
        }
    }

    #[test]
    fn polarization_skew_wraps_quarter_turn_references() {
        // A reference rotated by ±90° lands at the ends of (-π/2, π/2].
        let f = LocalFrame::new(&st(40.0, 0.0));
        let sat = geo_satellite_position(0.0);
        let k = (sat - f.origin).normalize();
        let horizontal = k.cross(f.up()).normalize();
        for r in [horizontal, -horizontal] {
            let a = polarization_skew(&f, sat, r).unwrap();
            assert_abs_diff_eq!(a.abs(), FRAC_PI_2, epsilon = 1e-9);
        }
    }

    #[test]
    fn gso_arc_separation_respects_elevation_mask() {
        let f = LocalFrame::new(&st(60.0, 0.0));
        let target = geo_satellite_position(0.0);
        let open = gso_arc_separation(&f, target, 0.0).unwrap();
        assert!(open.angle_rad < 1e-6);
        // A mask above the highest point of the arc hides it entirely.
        assert!(gso_arc_separation(&f, target, 30f64.to_radians()).is_none());
    }

    #[test]
    fn golden_section_finds_parabola_minimum() {
        let x = golden_section_min(-1.0, 3.0, &|x: f64| (x - 1.25).powi(2));
        assert_abs_diff_eq!(x, 1.25, epsilon = 1e-9);
    }

    #[test]
    fn coarse_scan_returns_none_when_all_slots_masked() {
        assert!(coarse_arc_scan(0.0, &|_| f64::INFINITY).is_none());
        let (s, lon) = coarse_arc_scan(0.0, &|l: f64| (l - 0.1).abs()).unwrap();
        assert!(s <= GSO_SCAN_STEP_RAD);
        assert_abs_diff_eq!(lon, 0.1, epsilon = GSO_SCAN_STEP_RAD);
    }

    #[test]
    fn horizon_mask_wraps_before_first_sample() {
        let m = HorizonMask::new(&[(0.5, 0.2), (PI, 0.0)]).unwrap();
        // 0.1 rad lies between the last sample (π) and the first one (0.5 + 2π).
        let t = (0.1 + TAU - PI) / (0.5 + TAU - PI);
        assert_abs_diff_eq!(m.min_elevation(0.1), 0.2 * t, epsilon = 1e-12);
    }
}
