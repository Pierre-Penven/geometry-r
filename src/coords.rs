//! Coordinate systems: geodetic (lat/lon/alt), Earth-Centred Earth-Fixed
//! (ECEF) and local topocentric East-North-Up (ENU) frames.

use crate::ellipsoid::{Ellipsoid, WGS84};
use crate::error::{GeoError, Result};
use crate::math::{Mat3, Vec3, wrap_pi, wrap_two_pi};

/// An ECEF position in metres. Alias of [`Vec3`] so that vector algebra
/// applies directly.
pub type Ecef = Vec3;

/// A geodetic position: latitude and longitude in radians on a reference
/// ellipsoid, and height above that ellipsoid in metres.
///
/// Longitude is always normalised to `(-π, π]`, so a point given at 190°E is
/// stored as -170°.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Geodetic {
    /// Geodetic latitude in radians, `[-π/2, π/2]`.
    pub lat_rad: f64,
    /// Longitude in radians, `(-π, π]`.
    pub lon_rad: f64,
    /// Height above the ellipsoid in metres.
    pub alt_m: f64,
}

impl Geodetic {
    /// Builds a geodetic position from radians.
    ///
    /// # Errors
    /// [`GeoError::InvalidCoordinate`] when a value is not finite or the
    /// latitude is outside `[-π/2, π/2]`.
    pub fn from_radians(lat_rad: f64, lon_rad: f64, alt_m: f64) -> Result<Self> {
        if !(lat_rad.is_finite() && lon_rad.is_finite() && alt_m.is_finite()) {
            return Err(GeoError::InvalidCoordinate("non-finite value".into()));
        }
        if lat_rad.abs() > std::f64::consts::FRAC_PI_2 + 1e-12 {
            return Err(GeoError::InvalidCoordinate(format!(
                "latitude {} deg outside [-90, 90]",
                lat_rad.to_degrees()
            )));
        }
        Ok(Self {
            lat_rad: lat_rad.clamp(-std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2),
            lon_rad: wrap_pi(lon_rad),
            alt_m,
        })
    }

    /// Builds a geodetic position from degrees (latitude, longitude) and
    /// metres (altitude).
    ///
    /// # Errors
    /// See [`Geodetic::from_radians`].
    pub fn from_degrees(lat_deg: f64, lon_deg: f64, alt_m: f64) -> Result<Self> {
        Self::from_radians(lat_deg.to_radians(), lon_deg.to_radians(), alt_m)
    }

    /// Latitude in degrees.
    #[must_use]
    pub fn lat_deg(&self) -> f64 {
        self.lat_rad.to_degrees()
    }

    /// Longitude in degrees, `(-180, 180]`.
    #[must_use]
    pub fn lon_deg(&self) -> f64 {
        self.lon_rad.to_degrees()
    }

    /// Converts to ECEF on WGS-84.
    #[must_use]
    pub fn to_ecef(&self) -> Ecef {
        geodetic_to_ecef(&WGS84, self)
    }

    /// Converts an ECEF position to geodetic on WGS-84.
    #[must_use]
    pub fn from_ecef(p: Ecef) -> Self {
        ecef_to_geodetic(&WGS84, p)
    }

    /// Unit vector of the geocentric direction of this point on a sphere
    /// (latitude interpreted as spherical). This is the representation used
    /// by the spherical-geometry functions.
    #[must_use]
    pub fn to_unit_vector(&self) -> Vec3 {
        Ellipsoid::normal(self.lat_rad, self.lon_rad)
    }
}

/// Converts geodetic coordinates to ECEF on the given ellipsoid.
#[inline]
#[must_use]
pub fn geodetic_to_ecef(ell: &Ellipsoid, g: &Geodetic) -> Ecef {
    let (sp, cp) = g.lat_rad.sin_cos();
    let (sl, cl) = g.lon_rad.sin_cos();
    let n = ell.a / (1.0 - ell.e2 * sp * sp).sqrt();
    let r = (n + g.alt_m) * cp;
    Vec3::new(r * cl, r * sl, (n * (1.0 - ell.e2) + g.alt_m) * sp)
}

/// Converts ECEF coordinates to geodetic on the given ellipsoid.
///
/// Uses Heikkinen's closed-form solution: no iteration, sub-millimetre
/// accuracy from the Earth's surface out to beyond the geostationary orbit.
#[must_use]
pub fn ecef_to_geodetic(ell: &Ellipsoid, p: Ecef) -> Geodetic {
    let lon = p.y.atan2(p.x);
    let rho2 = p.x * p.x + p.y * p.y;
    let rho = rho2.sqrt();

    if ell.e2 == 0.0 {
        let r = p.norm();
        return Geodetic {
            lat_rad: p.z.atan2(rho),
            lon_rad: lon,
            alt_m: r - ell.a,
        };
    }

    let (a, b, e2) = (ell.a, ell.b, ell.e2);
    let (a2, b2) = (a * a, b * b);
    let z2 = p.z * p.z;
    let e4 = e2 * e2;
    let f = 54.0 * b2 * z2;
    let g = rho2 + (1.0 - e2) * z2 - e2 * (a2 - b2);
    let c = e4 * f * rho2 / (g * g * g);
    let s = (1.0 + c + (c * c + 2.0 * c).sqrt()).cbrt();
    let k = s + 1.0 + 1.0 / s;
    let pp = f / (3.0 * k * k * g * g);
    let q = (1.0 + 2.0 * e4 * pp).sqrt();
    let r0_arg = 0.5 * a2 * (1.0 + 1.0 / q) - pp * (1.0 - e2) * z2 / (q * (1.0 + q)) - 0.5 * pp * rho2;
    let r0 = -(pp * e2 * rho) / (1.0 + q) + r0_arg.max(0.0).sqrt();
    let t = rho - e2 * r0;
    let u = (t * t + z2).sqrt();
    let v = (t * t + (1.0 - e2) * z2).sqrt();
    let z0 = b2 * p.z / (a * v);
    let alt = u * (1.0 - b2 / (a * v));
    let lat = (p.z + ell.ep2 * z0).atan2(rho);
    Geodetic {
        lat_rad: lat,
        lon_rad: lon,
        alt_m: alt,
    }
}

/// Azimuth / elevation / range of a target seen from an observer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LookAngles {
    /// Azimuth in radians, clockwise from geodetic north, `[0, 2π)`.
    pub azimuth_rad: f64,
    /// Elevation above the local horizontal plane (ellipsoid tangent) in
    /// radians, `[-π/2, π/2]`.
    pub elevation_rad: f64,
    /// Slant range in metres.
    pub range_m: f64,
}

/// A local topocentric East-North-Up frame anchored at a point.
///
/// Building the frame costs a few trigonometric calls; transforming points
/// afterwards is a 3×3 product. Build it once per ground station and reuse it
/// for every satellite or epoch.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LocalFrame {
    /// Origin of the frame in ECEF.
    pub origin: Ecef,
    /// Geodetic coordinates of the origin.
    pub geodetic: Geodetic,
    /// Rotation ECEF → ENU (rows are the east, north and up unit vectors).
    pub rot: Mat3,
}

impl LocalFrame {
    /// Creates an ENU frame at a geodetic position on WGS-84.
    #[must_use]
    pub fn new(origin: &Geodetic) -> Self {
        Self::with_ellipsoid(&WGS84, origin)
    }

    /// Creates an ENU frame at a geodetic position on the given ellipsoid.
    #[must_use]
    pub fn with_ellipsoid(ell: &Ellipsoid, origin: &Geodetic) -> Self {
        let (sp, cp) = origin.lat_rad.sin_cos();
        let (sl, cl) = origin.lon_rad.sin_cos();
        let east = Vec3::new(-sl, cl, 0.0);
        let north = Vec3::new(-sp * cl, -sp * sl, cp);
        let up = Vec3::new(cp * cl, cp * sl, sp);
        Self {
            origin: geodetic_to_ecef(ell, origin),
            geodetic: *origin,
            rot: Mat3::from_rows(east, north, up),
        }
    }

    /// East unit vector in ECEF.
    #[must_use]
    pub fn east(&self) -> Vec3 {
        self.rot.rows[0]
    }

    /// North unit vector in ECEF.
    #[must_use]
    pub fn north(&self) -> Vec3 {
        self.rot.rows[1]
    }

    /// Up unit vector (ellipsoid normal) in ECEF.
    #[must_use]
    pub fn up(&self) -> Vec3 {
        self.rot.rows[2]
    }

    /// ECEF position → ENU coordinates (metres) relative to the origin.
    #[inline]
    #[must_use]
    pub fn to_enu(&self, p: Ecef) -> Vec3 {
        self.rot.mul_vec(p - self.origin)
    }

    /// ENU coordinates (metres) → ECEF position.
    #[inline]
    #[must_use]
    pub fn from_enu(&self, enu: Vec3) -> Ecef {
        self.rot.mul_vec_t(enu) + self.origin
    }

    /// ECEF direction (e.g. a velocity) → ENU direction. No translation.
    #[inline]
    #[must_use]
    pub fn dir_to_enu(&self, d: Vec3) -> Vec3 {
        self.rot.mul_vec(d)
    }

    /// Azimuth, elevation and range of an ECEF target from the origin.
    #[inline]
    #[must_use]
    pub fn look_angles(&self, target: Ecef) -> LookAngles {
        enu_to_look(self.to_enu(target))
    }

    /// ECEF unit direction of a pointing given by azimuth and elevation.
    #[must_use]
    pub fn direction(&self, azimuth_rad: f64, elevation_rad: f64) -> Vec3 {
        let (sa, ca) = azimuth_rad.sin_cos();
        let (se, ce) = elevation_rad.sin_cos();
        self.rot.mul_vec_t(Vec3::new(ce * sa, ce * ca, se))
    }
}

/// Converts an ENU vector to azimuth/elevation/range.
#[inline]
#[must_use]
pub fn enu_to_look(enu: Vec3) -> LookAngles {
    let horiz = enu.x.hypot(enu.y);
    LookAngles {
        azimuth_rad: wrap_two_pi(enu.x.atan2(enu.y)),
        elevation_rad: enu.z.atan2(horiz),
        range_m: enu.norm(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn geodetic_validation() {
        assert!(matches!(
            Geodetic::from_degrees(91.0, 0.0, 0.0),
            Err(GeoError::InvalidCoordinate(_))
        ));
        assert!(matches!(
            Geodetic::from_degrees(f64::NAN, 0.0, 0.0),
            Err(GeoError::InvalidCoordinate(_))
        ));
        assert!(matches!(
            Geodetic::from_degrees(0.0, 0.0, f64::INFINITY),
            Err(GeoError::InvalidCoordinate(_))
        ));
        let g = Geodetic::from_degrees(10.0, 190.0, 5.0).unwrap();
        assert_abs_diff_eq!(g.lon_deg(), -170.0, epsilon = 1e-12);
        assert_abs_diff_eq!(g.lat_deg(), 10.0, epsilon = 1e-12);
    }

    #[test]
    fn known_ecef_values() {
        let p = Geodetic::from_degrees(0.0, 0.0, 0.0).unwrap().to_ecef();
        assert_abs_diff_eq!(p.x, WGS84.a, epsilon = 1e-6);
        let p = Geodetic::from_degrees(90.0, 0.0, 0.0).unwrap().to_ecef();
        assert_abs_diff_eq!(p.z, WGS84.b, epsilon = 1e-6);
        assert_abs_diff_eq!(p.x, 0.0, epsilon = 1e-6);
    }

    #[test]
    fn round_trip_surface_to_geo() {
        for &lat in &[-90.0, -89.999, -45.0, -1e-9, 0.0, 12.3, 60.0, 89.9999, 90.0] {
            for &lon in &[-180.0, -179.9, -30.0, 0.0, 77.7, 179.99, 180.0] {
                for &alt in &[-400.0, 0.0, 1500.0, 550e3, 35_786e3, 400_000e3] {
                    let g = Geodetic::from_degrees(lat, lon, alt).unwrap();
                    let back = Geodetic::from_ecef(g.to_ecef());
                    assert_abs_diff_eq!(back.lat_rad, g.lat_rad, epsilon = 1e-11);
                    assert_abs_diff_eq!(back.alt_m, g.alt_m, epsilon = 1e-3);
                    if lat.abs() < 90.0 {
                        assert_abs_diff_eq!(crate::math::wrap_pi(back.lon_rad - g.lon_rad), 0.0, epsilon = 1e-12);
                    }
                }
            }
        }
    }

    #[test]
    fn sphere_round_trip() {
        let s = Ellipsoid::sphere(6_371_000.0);
        let g = Geodetic::from_degrees(33.0, -100.0, 1234.0).unwrap();
        let back = ecef_to_geodetic(&s, geodetic_to_ecef(&s, &g));
        assert_abs_diff_eq!(back.lat_rad, g.lat_rad, epsilon = 1e-12);
        assert_abs_diff_eq!(back.alt_m, 1234.0, epsilon = 1e-6);
    }

    #[test]
    fn enu_and_look_angles() {
        let st = Geodetic::from_degrees(45.0, 7.0, 300.0).unwrap();
        let f = LocalFrame::new(&st);
        // Point straight up.
        let up = f.from_enu(Vec3::new(0.0, 0.0, 1000.0));
        let la = f.look_angles(up);
        assert_abs_diff_eq!(la.elevation_rad, std::f64::consts::FRAC_PI_2, epsilon = 1e-9);
        assert_abs_diff_eq!(la.range_m, 1000.0, epsilon = 1e-6);
        // Point due east on the horizon.
        let la = f.look_angles(f.from_enu(Vec3::new(10.0, 0.0, 0.0)));
        assert_abs_diff_eq!(la.azimuth_rad, std::f64::consts::FRAC_PI_2, epsilon = 1e-9);
        assert_abs_diff_eq!(la.elevation_rad, 0.0, epsilon = 1e-9);
        // West gives 270°.
        let la = f.look_angles(f.from_enu(Vec3::new(-10.0, 0.0, 0.0)));
        assert_abs_diff_eq!(la.azimuth_rad.to_degrees(), 270.0, epsilon = 1e-9);
        // Direction is consistent with look angles.
        let d = f.direction(0.3, 0.2);
        let la = f.look_angles(f.origin + d * 1e6);
        assert_abs_diff_eq!(la.azimuth_rad, 0.3, epsilon = 1e-9);
        assert_abs_diff_eq!(la.elevation_rad, 0.2, epsilon = 1e-9);
        assert_abs_diff_eq!(f.east().dot(f.north()), 0.0, epsilon = 1e-15);
        assert_abs_diff_eq!(f.up().dot(st.to_unit_vector()), 1.0, epsilon = 1e-15);
        assert_abs_diff_eq!(f.dir_to_enu(f.up()).z, 1.0, epsilon = 1e-15);
    }
}
