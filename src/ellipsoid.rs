//! Reference ellipsoids and the geometric primitives that depend on them
//! (radii of curvature, surface normals, ray intersection).

use crate::math::Vec3;

/// An oblate ellipsoid of revolution about the Z axis.
///
/// Derived quantities are pre-computed so that hot paths never recompute
/// them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ellipsoid {
    /// Semi-major (equatorial) axis in metres.
    pub a: f64,
    /// Flattening `(a − b) / a`.
    pub f: f64,
    /// Semi-minor (polar) axis in metres.
    pub b: f64,
    /// First eccentricity squared `(a² − b²) / a²`.
    pub e2: f64,
    /// Second eccentricity squared `(a² − b²) / b²`.
    pub ep2: f64,
}

/// WGS-84 ellipsoid (GPS, GeoJSON and most satellite systems).
pub const WGS84: Ellipsoid = Ellipsoid {
    a: 6_378_137.0,
    f: 1.0 / 298.257_223_563,
    b: 6_356_752.314_245_179,
    e2: 0.006_694_379_990_141_317,
    ep2: 0.006_739_496_742_276_435,
};

/// Mean Earth radius `R1 = (2a + b) / 3` for WGS-84, in metres.
pub const EARTH_MEAN_RADIUS_M: f64 = 6_371_008.771_415;

/// Authalic (equal-area) radius of WGS-84, in metres. Spherical areas
/// computed with this radius match ellipsoidal areas of large regions closely.
pub const EARTH_AUTHALIC_RADIUS_M: f64 = 6_371_007.180_918;

/// Earth rotation rate (IERS), in rad/s.
pub const EARTH_ROTATION_RATE_RAD_S: f64 = 7.292_115_146_706_979e-5;

/// Standard gravitational parameter of the Earth (WGS-84), in m³/s².
pub const EARTH_GM_M3_S2: f64 = 3.986_004_418e14;

/// Radius of the geostationary orbit, in metres.
pub const GEO_RADIUS_M: f64 = 42_164_172.0;

/// Speed of light in vacuum, in m/s.
pub const SPEED_OF_LIGHT_M_S: f64 = 299_792_458.0;

impl Ellipsoid {
    /// Builds an ellipsoid from its semi-major axis (m) and flattening.
    #[must_use]
    pub fn new(a: f64, f: f64) -> Self {
        let b = a * (1.0 - f);
        let e2 = f * (2.0 - f);
        let ep2 = e2 / ((1.0 - f) * (1.0 - f));
        Self { a, f, b, e2, ep2 }
    }

    /// A perfect sphere of radius `r` metres.
    #[must_use]
    pub const fn sphere(r: f64) -> Self {
        Self {
            a: r,
            f: 0.0,
            b: r,
            e2: 0.0,
            ep2: 0.0,
        }
    }

    /// True when the ellipsoid is a sphere (zero flattening).
    #[must_use]
    pub fn is_sphere(&self) -> bool {
        self.f == 0.0
    }

    /// Prime-vertical radius of curvature `N(φ)` in metres.
    #[inline]
    #[must_use]
    pub fn prime_vertical_radius(&self, lat_rad: f64) -> f64 {
        let s = lat_rad.sin();
        self.a / (1.0 - self.e2 * s * s).sqrt()
    }

    /// Meridional radius of curvature `M(φ)` in metres.
    #[inline]
    #[must_use]
    pub fn meridional_radius(&self, lat_rad: f64) -> f64 {
        let s = lat_rad.sin();
        let w2 = 1.0 - self.e2 * s * s;
        self.a * (1.0 - self.e2) / (w2 * w2.sqrt())
    }

    /// Distance from the Earth's centre to the surface at geodetic latitude
    /// `lat_rad`, in metres.
    #[must_use]
    pub fn geocentric_radius(&self, lat_rad: f64) -> f64 {
        let (s, c) = lat_rad.sin_cos();
        let (a2, b2) = (self.a * self.a, self.b * self.b);
        let num = (a2 * c).powi(2) + (b2 * s).powi(2);
        let den = (self.a * c).powi(2) + (self.b * s).powi(2);
        (num / den).sqrt()
    }

    /// Outward unit normal of the ellipsoid at the given geodetic latitude
    /// and longitude (this is the local "up" direction).
    #[inline]
    #[must_use]
    pub fn normal(lat_rad: f64, lon_rad: f64) -> Vec3 {
        let (sp, cp) = lat_rad.sin_cos();
        let (sl, cl) = lon_rad.sin_cos();
        Vec3::new(cp * cl, cp * sl, sp)
    }

    /// Distance `t ≥ 0` along the ray `origin + t·dir` to the first
    /// intersection with the ellipsoid surface (raised by `height_m`), or
    /// `None` if the ray misses it or points away from it.
    ///
    /// `dir` need not be normalised; `t` is expressed in units of `|dir|`.
    /// When `origin` is inside the surface, the exit point is returned.
    #[must_use]
    pub fn ray_intersection(&self, origin: Vec3, dir: Vec3, height_m: f64) -> Option<f64> {
        let (a, b) = (self.a + height_m, self.b + height_m);
        // Scale space so the ellipsoid becomes the unit sphere.
        let o = Vec3::new(origin.x / a, origin.y / a, origin.z / b);
        let d = Vec3::new(dir.x / a, dir.y / a, dir.z / b);
        let qa = d.norm2();
        if qa == 0.0 {
            return None;
        }
        let qb = o.dot(d);
        let qc = o.norm2() - 1.0;
        let disc = qb * qb - qa * qc;
        if disc < 0.0 {
            return None;
        }
        let sq = disc.sqrt();
        // Numerically stable roots.
        let q = -(qb + qb.signum() * sq);
        let (mut t0, mut t1) = if q != 0.0 { (q / qa, qc / q) } else { (0.0, 0.0) };
        if t0 > t1 {
            std::mem::swap(&mut t0, &mut t1);
        }
        if t0 >= 0.0 {
            Some(t0)
        } else if t1 >= 0.0 {
            Some(t1)
        } else {
            None
        }
    }

    /// True when the straight segment between two ECEF points passes below
    /// the ellipsoid raised by `grazing_height_m` (i.e. the Earth blocks the
    /// line of sight). Endpoints themselves may be on or below the surface:
    /// only the open interior of the segment is tested.
    #[must_use]
    pub fn segment_blocked(&self, p: Vec3, q: Vec3, grazing_height_m: f64) -> bool {
        let (a, b) = (self.a + grazing_height_m, self.b + grazing_height_m);
        let ps = Vec3::new(p.x / a, p.y / a, p.z / b);
        let qs = Vec3::new(q.x / a, q.y / a, q.z / b);
        let d = qs - ps;
        let dd = d.norm2();
        if dd == 0.0 {
            return false;
        }
        // Closest approach of the scaled segment to the origin.
        let t = (-ps.dot(d) / dd).clamp(0.0, 1.0);
        if t <= 0.0 || t >= 1.0 {
            return false;
        }
        // Small tolerance so that a point exactly on the surface looking up
        // is not considered blocked by itself.
        (ps + d * t).norm2() < 1.0 - 1e-12
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn wgs84_matches_derived_values() {
        let d = Ellipsoid::new(WGS84.a, WGS84.f);
        assert_abs_diff_eq!(d.b, WGS84.b, epsilon = 1e-8);
        assert_abs_diff_eq!(d.e2, WGS84.e2, epsilon = 1e-17);
        assert_abs_diff_eq!(d.ep2, WGS84.ep2, epsilon = 1e-17);
    }

    #[test]
    fn wgs84_constants() {
        assert_abs_diff_eq!(WGS84.b, 6_356_752.314_245, epsilon = 1e-5);
        assert_abs_diff_eq!(WGS84.e2, 0.006_694_379_990_14, epsilon = 1e-14);
        assert_abs_diff_eq!(WGS84.ep2, 0.006_739_496_742_28, epsilon = 1e-14);
        assert!(!WGS84.is_sphere());
        assert!(Ellipsoid::sphere(1.0).is_sphere());
    }

    #[test]
    fn radii() {
        assert_abs_diff_eq!(WGS84.geocentric_radius(0.0), WGS84.a, epsilon = 1e-6);
        assert_abs_diff_eq!(
            WGS84.geocentric_radius(std::f64::consts::FRAC_PI_2),
            WGS84.b,
            epsilon = 1e-6
        );
        assert_abs_diff_eq!(WGS84.prime_vertical_radius(0.0), WGS84.a, epsilon = 1e-6);
        assert_abs_diff_eq!(WGS84.meridional_radius(0.0), WGS84.a * (1.0 - WGS84.e2), epsilon = 1e-6);
    }

    #[test]
    fn ray_hits_and_misses() {
        let o = Vec3::new(10.0 * WGS84.a, 0.0, 0.0);
        let t = WGS84.ray_intersection(o, -Vec3::X, 0.0).expect("hits");
        assert_abs_diff_eq!(t, 9.0 * WGS84.a, epsilon = 1e-6);
        assert!(WGS84.ray_intersection(o, Vec3::X, 0.0).is_none());
        assert!(WGS84.ray_intersection(o, Vec3::Y, 0.0).is_none());
        assert!(WGS84.ray_intersection(o, Vec3::ZERO, 0.0).is_none());
        // From the centre the exit point is found.
        let t = WGS84.ray_intersection(Vec3::ZERO, Vec3::Z, 0.0).expect("exit");
        assert_abs_diff_eq!(t, WGS84.b, epsilon = 1e-6);
        // Raised surface.
        let t = WGS84.ray_intersection(o, -Vec3::X, 1000.0).expect("hits");
        assert_abs_diff_eq!(t, 9.0 * WGS84.a - 1000.0, epsilon = 1e-6);
    }

    #[test]
    fn segment_blocking() {
        let a = Vec3::new(WGS84.a + 500e3, 0.0, 0.0);
        let b = Vec3::new(-(WGS84.a + 500e3), 0.0, 0.0);
        assert!(WGS84.segment_blocked(a, b, 0.0));
        let c = Vec3::new(WGS84.a + 500e3, 1000e3, 0.0);
        assert!(!WGS84.segment_blocked(a, c, 0.0));
        assert!(!WGS84.segment_blocked(a, a, 0.0));
        // A ground point looking straight up is not blocked by itself.
        let g = Vec3::new(WGS84.a, 0.0, 0.0);
        assert!(!WGS84.segment_blocked(g, a, 0.0));
    }

    #[test]
    fn tangent_ray_from_surface_point_hits_at_origin() {
        let o = Vec3::new(WGS84.a, 0.0, 0.0);
        assert_eq!(WGS84.ray_intersection(o, Vec3::Y, 0.0), Some(0.0));
    }
}
