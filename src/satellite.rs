//! Satellite-centred geometry: body frames, scan angles, inverse scan to the
//! ground, and the closed-form nadir / elevation / central-angle relations.

use crate::coords::{Ecef, Geodetic, ecef_to_geodetic};
use crate::ellipsoid::{EARTH_ROTATION_RATE_RAD_S, Ellipsoid, GEO_RADIUS_M, WGS84};
use crate::error::{GeoError, Result};
use crate::math::{Mat3, Vec3, wrap_two_pi};

/// Which "down" direction defines the nadir of a satellite frame.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum NadirReference {
    /// Towards the Earth's centre.
    #[default]
    Geocentric,
    /// Along the ellipsoid normal through the sub-satellite point.
    Geodetic,
}

/// Orthonormal body frame of a satellite (or any platform).
///
/// Body axes are expressed in ECEF. For the frames built here, `+Z` points
/// to nadir (antenna boresight of a nadir-pointing payload), `+X` roughly
/// along the velocity, and `+Y` completes the right-handed triad (towards
/// the negative orbit normal), i.e. the classic LVLH / roll-pitch-yaw axes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SatelliteFrame {
    /// Platform position in ECEF (metres).
    pub position: Ecef,
    /// Rotation ECEF → body. Rows are the body X, Y, Z axes in ECEF.
    pub rot: Mat3,
}

impl SatelliteFrame {
    /// Builds a frame from explicit body axes given in ECEF. The axes are
    /// re-orthonormalised (Z kept, X made orthogonal to it).
    ///
    /// # Errors
    /// [`GeoError::InvalidArgument`] if the axes are degenerate or parallel.
    pub fn from_axes(position: Ecef, x_axis: Vec3, z_axis: Vec3) -> Result<Self> {
        let z = z_axis
            .try_normalize()
            .ok_or_else(|| GeoError::InvalidArgument("zero z axis".into()))?;
        let x = (x_axis - z * x_axis.dot(z))
            .try_normalize()
            .filter(|_| x_axis.cross(z).norm() > 1e-12 * x_axis.norm())
            .ok_or_else(|| GeoError::InvalidArgument("x axis parallel to z axis".into()))?;
        let y = z.cross(x);
        Ok(Self {
            position,
            rot: Mat3::from_rows(x, y, z),
        })
    }

    /// Local-vertical local-horizontal frame from position and velocity
    /// (both ECEF; use the inertial or the Earth-fixed velocity consistently
    /// with how the attitude is defined).
    ///
    /// `Z = nadir`, `Y = −(r × v)` (negative orbit normal), `X = Y × Z`.
    ///
    /// # Errors
    /// [`GeoError::InvalidArgument`] if the position is zero or the velocity
    /// is parallel to it.
    pub fn lvlh(position: Ecef, velocity: Vec3, nadir: NadirReference) -> Result<Self> {
        let z = nadir_direction(position, nadir)?;
        let h = position.cross(velocity);
        let y = (-h)
            .try_normalize()
            .ok_or_else(|| GeoError::InvalidArgument("velocity parallel to position".into()))?;
        // Re-orthogonalise Y against Z. Y is exactly orthogonal to the radial
        // direction and the geodetic nadir is within 0.2° of it, so the
        // projection never vanishes.
        let y = (y - z * y.dot(z)).normalize();
        let x = y.cross(z);
        Ok(Self {
            position,
            rot: Mat3::from_rows(x, y, z),
        })
    }

    /// Nadir-pointing frame of a geostationary satellite: the velocity is
    /// taken as the Earth-rotation velocity `ω × r`, so `X` points east, `Y`
    /// south and `Z` to nadir.
    ///
    /// # Errors
    /// [`GeoError::InvalidArgument`] if the position is on the polar axis.
    pub fn geo_nadir(position: Ecef, nadir: NadirReference) -> Result<Self> {
        let v = Vec3::new(
            -EARTH_ROTATION_RATE_RAD_S * position.y,
            EARTH_ROTATION_RATE_RAD_S * position.x,
            0.0,
        );
        Self::lvlh(position, v, nadir)
    }

    /// Applies an attitude bias: body axes rotated by `roll` about X, then
    /// `pitch` about Y, then `yaw` about Z (intrinsic Z-Y-X order, radians).
    #[must_use]
    pub fn with_attitude(&self, roll_rad: f64, pitch_rad: f64, yaw_rad: f64) -> Self {
        let r = Mat3::rot_z(yaw_rad)
            .mul_mat(&Mat3::rot_y(pitch_rad))
            .mul_mat(&Mat3::rot_x(roll_rad));
        Self {
            position: self.position,
            rot: r.transpose().mul_mat(&self.rot),
        }
    }

    /// Body axis X in ECEF.
    #[must_use]
    pub fn x_axis(&self) -> Vec3 {
        self.rot.rows[0]
    }

    /// Body axis Y in ECEF.
    #[must_use]
    pub fn y_axis(&self) -> Vec3 {
        self.rot.rows[1]
    }

    /// Body axis Z (boresight / nadir) in ECEF.
    #[must_use]
    pub fn z_axis(&self) -> Vec3 {
        self.rot.rows[2]
    }

    /// ECEF direction → body direction.
    #[inline]
    #[must_use]
    pub fn to_body(&self, dir_ecef: Vec3) -> Vec3 {
        self.rot.mul_vec(dir_ecef)
    }

    /// Body direction → ECEF direction.
    #[inline]
    #[must_use]
    pub fn from_body(&self, dir_body: Vec3) -> Vec3 {
        self.rot.mul_vec_t(dir_body)
    }

    /// Unit body-frame direction from the platform to an ECEF target.
    ///
    /// # Errors
    /// [`GeoError::InvalidArgument`] if the target coincides with the platform.
    pub fn direction_to(&self, target: Ecef) -> Result<Vec3> {
        self.to_body(target - self.position)
            .try_normalize()
            .ok_or_else(|| GeoError::InvalidArgument("target coincides with platform".into()))
    }

    /// Scan angles of an ECEF target in this body frame.
    ///
    /// # Errors
    /// See [`SatelliteFrame::direction_to`].
    pub fn scan_angles(&self, target: Ecef) -> Result<ScanAngles> {
        Ok(ScanAngles::from_direction(self.direction_to(target)?))
    }

    /// Ground point (on the ellipsoid raised by `height_m`) hit by the ray
    /// leaving the platform along `dir_body`.
    ///
    /// # Errors
    /// [`GeoError::NoIntersection`] when the ray misses the Earth.
    pub fn ground_intersection(&self, dir_body: Vec3, ell: &Ellipsoid, height_m: f64) -> Result<Ecef> {
        let d = self.from_body(dir_body);
        ell.ray_intersection(self.position, d, height_m)
            .map(|t| self.position + d * t)
            .ok_or(GeoError::NoIntersection)
    }

    /// Ground point on WGS-84 for scan angles.
    ///
    /// # Errors
    /// [`GeoError::NoIntersection`] when the ray misses the Earth.
    pub fn scan_to_ground(&self, scan: &ScanAngles) -> Result<Geodetic> {
        let p = self.ground_intersection(scan.direction(), &WGS84, 0.0)?;
        Ok(ecef_to_geodetic(&WGS84, p))
    }
}

/// A body-frame pointing direction expressed in the usual antenna
/// conventions. All three are views of the same unit vector `d`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScanAngles {
    /// Off-boresight angle θ = angle between `d` and `+Z`, `[0, π]`.
    pub theta_rad: f64,
    /// Azimuth φ of `d` in the X-Y plane, from `+X` towards `+Y`, `[0, 2π)`.
    pub phi_rad: f64,
    /// Direction cosine u = d·X.
    pub u: f64,
    /// Direction cosine v = d·Y.
    pub v: f64,
    /// Tilt about −Y seen in the X-Z plane: `atan2(dx, dz)`.
    pub tilt_x_rad: f64,
    /// Tilt about X seen in the Y-Z plane: `atan2(dy, dz)`.
    pub tilt_y_rad: f64,
}

impl ScanAngles {
    /// Computes every convention from a body direction (need not be unit).
    #[must_use]
    pub fn from_direction(d: Vec3) -> Self {
        let d = d.normalize();
        Self {
            theta_rad: d.x.hypot(d.y).atan2(d.z),
            phi_rad: wrap_two_pi(d.y.atan2(d.x)),
            u: d.x,
            v: d.y,
            tilt_x_rad: d.x.atan2(d.z),
            tilt_y_rad: d.y.atan2(d.z),
        }
    }

    /// From θ (off-boresight) and φ (azimuth from +X towards +Y).
    #[must_use]
    pub fn from_theta_phi(theta_rad: f64, phi_rad: f64) -> Self {
        let (st, ct) = theta_rad.sin_cos();
        let (sp, cp) = phi_rad.sin_cos();
        Self::from_direction(Vec3::new(st * cp, st * sp, ct))
    }

    /// From direction cosines (forward hemisphere, `w = +√(1−u²−v²)`).
    ///
    /// # Errors
    /// [`GeoError::InvalidArgument`] when `u² + v² > 1`.
    pub fn from_uv(u: f64, v: f64) -> Result<Self> {
        let w2 = 1.0 - u * u - v * v;
        if w2.is_nan() || w2 < -1e-15 {
            return Err(GeoError::InvalidArgument(format!("u²+v² = {} > 1", 1.0 - w2)));
        }
        Ok(Self::from_direction(Vec3::new(u, v, w2.max(0.0).sqrt())))
    }

    /// From X/Y tilt angles (both in `(-π/2, π/2)`, forward hemisphere).
    ///
    /// # Errors
    /// [`GeoError::InvalidArgument`] when a tilt is outside `(-π/2, π/2)`.
    pub fn from_tilts(tilt_x_rad: f64, tilt_y_rad: f64) -> Result<Self> {
        let lim = std::f64::consts::FRAC_PI_2;
        if !(tilt_x_rad.abs() < lim && tilt_y_rad.abs() < lim) {
            return Err(GeoError::InvalidArgument("tilt outside (-90°, 90°)".into()));
        }
        Ok(Self::from_direction(Vec3::new(tilt_x_rad.tan(), tilt_y_rad.tan(), 1.0)))
    }

    /// Unit body direction.
    #[must_use]
    pub fn direction(&self) -> Vec3 {
        let (st, ct) = self.theta_rad.sin_cos();
        let (sp, cp) = self.phi_rad.sin_cos();
        Vec3::new(st * cp, st * sp, ct)
    }
}

/// Unit nadir direction of a platform.
///
/// # Errors
/// [`GeoError::InvalidArgument`] for a zero position.
pub fn nadir_direction(position: Ecef, nadir: NadirReference) -> Result<Vec3> {
    let n = match nadir {
        NadirReference::Geocentric => -position,
        NadirReference::Geodetic => {
            let g = ecef_to_geodetic(&WGS84, position);
            -Ellipsoid::normal(g.lat_rad, g.lon_rad)
        }
    };
    n.try_normalize()
        .ok_or_else(|| GeoError::InvalidArgument("zero position".into()))
}

/// ECEF position of a geostationary satellite at the given longitude
/// (radians), on the ideal equatorial orbit of radius [`GEO_RADIUS_M`].
#[must_use]
pub fn geo_satellite_position(lon_rad: f64) -> Ecef {
    let (s, c) = lon_rad.sin_cos();
    Vec3::new(GEO_RADIUS_M * c, GEO_RADIUS_M * s, 0.0)
}

/// Geodetic sub-satellite point (geodetic nadir, altitude zero on WGS-84).
#[must_use]
pub fn sub_satellite_point(position: Ecef) -> Geodetic {
    let mut g = ecef_to_geodetic(&WGS84, position);
    g.alt_m = 0.0;
    g
}

/// Off-nadir angle (radians) of a target seen from the satellite, relative
/// to the geocentric nadir.
#[must_use]
pub fn off_nadir_angle(sat: Ecef, target: Ecef) -> f64 {
    (-sat).angle_to(target - sat)
}

/// Earth central angle (radians) between the satellite and a target.
#[must_use]
pub fn earth_central_angle(sat: Ecef, target: Ecef) -> f64 {
    sat.angle_to(target)
}

/// Spherical-Earth nadir angle η seen from radius `sat_radius_m` for a
/// ground elevation ε: `sin η = R cos ε / r`.
#[must_use]
pub fn nadir_angle_for_elevation(sat_radius_m: f64, earth_radius_m: f64, elevation_rad: f64) -> f64 {
    (earth_radius_m * elevation_rad.cos() / sat_radius_m)
        .clamp(-1.0, 1.0)
        .asin()
}

/// Spherical-Earth ground elevation ε for a nadir angle η:
/// `cos ε = r sin η / R`. Returns `None` past the horizon.
#[must_use]
pub fn elevation_for_nadir_angle(sat_radius_m: f64, earth_radius_m: f64, nadir_rad: f64) -> Option<f64> {
    let c = sat_radius_m * nadir_rad.sin() / earth_radius_m;
    (c <= 1.0).then(|| c.max(-1.0).acos())
}

/// Spherical-Earth central angle λ = 90° − ε − η for a ground elevation.
#[must_use]
pub fn central_angle_for_elevation(sat_radius_m: f64, earth_radius_m: f64, elevation_rad: f64) -> f64 {
    std::f64::consts::FRAC_PI_2 - elevation_rad - nadir_angle_for_elevation(sat_radius_m, earth_radius_m, elevation_rad)
}

/// Spherical-Earth slant range at a ground elevation (law of cosines).
#[must_use]
pub fn slant_range_for_elevation(sat_radius_m: f64, earth_radius_m: f64, elevation_rad: f64) -> f64 {
    let lambda = central_angle_for_elevation(sat_radius_m, earth_radius_m, elevation_rad);
    (sat_radius_m * sat_radius_m + earth_radius_m * earth_radius_m - 2.0 * sat_radius_m * earth_radius_m * lambda.cos())
        .sqrt()
}

/// Nadir angle of the Earth's limb from radius `sat_radius_m`: `asin(R/r)`.
#[must_use]
pub fn horizon_nadir_angle(sat_radius_m: f64, earth_radius_m: f64) -> f64 {
    (earth_radius_m / sat_radius_m).clamp(-1.0, 1.0).asin()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coords::LocalFrame;
    use crate::ellipsoid::EARTH_MEAN_RADIUS_M;
    use approx::assert_abs_diff_eq;

    #[test]
    fn geo_frame_axes() {
        let f = SatelliteFrame::geo_nadir(geo_satellite_position(0.0), NadirReference::Geocentric).unwrap();
        assert_abs_diff_eq!(f.z_axis().x, -1.0, epsilon = 1e-12); // nadir
        assert_abs_diff_eq!(f.x_axis().y, 1.0, epsilon = 1e-12); // east
        assert_abs_diff_eq!(f.y_axis().z, -1.0, epsilon = 1e-12); // south
    }

    #[test]
    fn frame_errors() {
        assert!(matches!(
            SatelliteFrame::lvlh(Vec3::ZERO, Vec3::X, NadirReference::Geocentric),
            Err(GeoError::InvalidArgument(_))
        ));
        assert!(matches!(
            SatelliteFrame::lvlh(Vec3::X * 7e6, Vec3::X, NadirReference::Geocentric),
            Err(GeoError::InvalidArgument(_))
        ));
        assert!(matches!(
            SatelliteFrame::from_axes(Vec3::X, Vec3::Z, Vec3::Z),
            Err(GeoError::InvalidArgument(_))
        ));
        assert!(matches!(
            SatelliteFrame::from_axes(Vec3::X, Vec3::Z, Vec3::ZERO),
            Err(GeoError::InvalidArgument(_))
        ));
        let f = SatelliteFrame::from_axes(Vec3::X, Vec3::Y, Vec3::Z).unwrap();
        assert!(matches!(f.direction_to(Vec3::X), Err(GeoError::InvalidArgument(_))));
    }

    #[test]
    fn scan_conventions_round_trip() {
        let s = ScanAngles::from_theta_phi(0.1, 2.0);
        let d = s.direction();
        let s2 = ScanAngles::from_uv(s.u, s.v).unwrap();
        assert_abs_diff_eq!(s2.theta_rad, 0.1, epsilon = 1e-12);
        assert_abs_diff_eq!(s2.phi_rad, 2.0, epsilon = 1e-12);
        let s3 = ScanAngles::from_tilts(s.tilt_x_rad, s.tilt_y_rad).unwrap();
        assert_abs_diff_eq!(s3.direction().dot(d), 1.0, epsilon = 1e-12);
        assert!(matches!(
            ScanAngles::from_uv(0.8, 0.8),
            Err(GeoError::InvalidArgument(_))
        ));
        assert!(matches!(
            ScanAngles::from_tilts(1.6, 0.0),
            Err(GeoError::InvalidArgument(_))
        ));
        assert!(matches!(
            ScanAngles::from_uv(f64::NAN, 0.0),
            Err(GeoError::InvalidArgument(_))
        ));
    }

    #[test]
    fn scan_to_ground_round_trip() {
        let f =
            SatelliteFrame::geo_nadir(geo_satellite_position(10f64.to_radians()), NadirReference::Geocentric).unwrap();
        let target = Geodetic::from_degrees(45.0, 2.0, 0.0).unwrap();
        let s = f.scan_angles(target.to_ecef()).unwrap();
        assert!(s.theta_rad < 8.7f64.to_radians());
        let g = f.scan_to_ground(&s).unwrap();
        assert_abs_diff_eq!(g.lat_rad, target.lat_rad, epsilon = 1e-9);
        assert_abs_diff_eq!(g.lon_rad, target.lon_rad, epsilon = 1e-9);
        // Northern target from a GEO frame has negative v (Y points south).
        assert!(s.v < 0.0);
        // Pointing to space misses.
        let miss = ScanAngles::from_theta_phi(0.3, 0.0);
        assert!(matches!(f.scan_to_ground(&miss), Err(GeoError::NoIntersection)));
    }

    #[test]
    fn attitude_bias() {
        let f = SatelliteFrame::geo_nadir(geo_satellite_position(0.0), NadirReference::Geocentric).unwrap();
        let rolled = f.with_attitude(0.01, 0.0, 0.0);
        // Rolling about X tilts boresight (Z) by 0.01 rad.
        assert_abs_diff_eq!(rolled.z_axis().angle_to(f.z_axis()), 0.01, epsilon = 1e-12);
        assert_abs_diff_eq!(rolled.x_axis().dot(f.x_axis()), 1.0, epsilon = 1e-12);
        let yawed = f.with_attitude(0.0, 0.0, 0.2);
        assert_abs_diff_eq!(yawed.z_axis().dot(f.z_axis()), 1.0, epsilon = 1e-12);
        assert_abs_diff_eq!(yawed.x_axis().angle_to(f.x_axis()), 0.2, epsilon = 1e-12);
        // to_body / from_body inverse.
        let d = Vec3::new(0.3, -0.2, 0.9).normalize();
        assert_abs_diff_eq!(rolled.from_body(rolled.to_body(d)).dot(d), 1.0, epsilon = 1e-12);
    }

    #[test]
    fn geodetic_nadir_differs_slightly() {
        let p = Geodetic::from_degrees(45.0, 0.0, 500e3).unwrap().to_ecef();
        let gc = nadir_direction(p, NadirReference::Geocentric).unwrap();
        let gd = nadir_direction(p, NadirReference::Geodetic).unwrap();
        let ang = gc.angle_to(gd).to_degrees();
        assert!(ang > 0.1 && ang < 0.2, "{ang}");
        let f = SatelliteFrame::lvlh(p, Vec3::new(0.0, 7500.0, 0.0), NadirReference::Geodetic).unwrap();
        assert_abs_diff_eq!(f.z_axis().dot(gd), 1.0, epsilon = 1e-12);
        assert_abs_diff_eq!(f.x_axis().dot(f.z_axis()), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn closed_form_relations() {
        let r = GEO_RADIUS_M;
        let re = EARTH_MEAN_RADIUS_M;
        let eta = nadir_angle_for_elevation(r, re, 0.0);
        assert_abs_diff_eq!(eta, horizon_nadir_angle(r, re), epsilon = 1e-12);
        assert_abs_diff_eq!(eta.to_degrees(), 8.69, epsilon = 0.01);
        let lambda = central_angle_for_elevation(r, re, 0.0);
        assert_abs_diff_eq!(lambda.to_degrees(), 81.3, epsilon = 0.05);
        let el = elevation_for_nadir_angle(r, re, nadir_angle_for_elevation(r, re, 0.3)).unwrap();
        assert_abs_diff_eq!(el, 0.3, epsilon = 1e-9);
        assert!(elevation_for_nadir_angle(r, re, 0.2).is_none());
        assert_abs_diff_eq!(
            slant_range_for_elevation(r, re, std::f64::consts::FRAC_PI_2),
            r - re,
            epsilon = 1e-3
        );
    }

    #[test]
    fn nadir_and_central_angles() {
        let sat = geo_satellite_position(0.0);
        let ssp = sub_satellite_point(sat);
        assert_abs_diff_eq!(ssp.lat_rad, 0.0, epsilon = 1e-12);
        assert_eq!(ssp.alt_m, 0.0);
        let t = Geodetic::from_degrees(0.0, 30.0, 0.0).unwrap().to_ecef();
        assert_abs_diff_eq!(earth_central_angle(sat, t).to_degrees(), 30.0, epsilon = 1e-9);
        let eta = off_nadir_angle(sat, t);
        // Consistent with the elevation seen from the ground.
        let el = LocalFrame::new(&Geodetic::from_degrees(0.0, 30.0, 0.0).unwrap())
            .look_angles(sat)
            .elevation_rad;
        let eta2 = nadir_angle_for_elevation(GEO_RADIUS_M, crate::ellipsoid::WGS84.a, el);
        assert_abs_diff_eq!(eta, eta2, epsilon = 1e-9);
    }
}
