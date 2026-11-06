//! Small, allocation-free 3D vector and matrix primitives.
//!
//! Everything in the crate is built on these two types. They are `Copy`,
//! `#[inline]` and free of bounds checks so that the compiler can vectorise
//! hot loops (batch look-angle or footprint computations).

use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

/// A 3D vector of `f64`. Used for ECEF positions (metres), velocities (m/s)
/// and unit direction vectors.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vec3 {
    /// X component.
    pub x: f64,
    /// Y component.
    pub y: f64,
    /// Z component.
    pub z: f64,
}

impl Vec3 {
    /// The zero vector.
    pub const ZERO: Self = Self::new(0.0, 0.0, 0.0);
    /// Unit X axis.
    pub const X: Self = Self::new(1.0, 0.0, 0.0);
    /// Unit Y axis.
    pub const Y: Self = Self::new(0.0, 1.0, 0.0);
    /// Unit Z axis.
    pub const Z: Self = Self::new(0.0, 0.0, 1.0);

    /// Creates a vector from its components.
    #[inline]
    #[must_use]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Dot product.
    #[inline]
    #[must_use]
    pub fn dot(self, o: Self) -> f64 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    /// Cross product.
    #[inline]
    #[must_use]
    pub fn cross(self, o: Self) -> Self {
        Self::new(
            self.y * o.z - self.z * o.y,
            self.z * o.x - self.x * o.z,
            self.x * o.y - self.y * o.x,
        )
    }

    /// Squared Euclidean norm.
    #[inline]
    #[must_use]
    pub fn norm2(self) -> f64 {
        self.dot(self)
    }

    /// Euclidean norm.
    #[inline]
    #[must_use]
    pub fn norm(self) -> f64 {
        self.norm2().sqrt()
    }

    /// Returns the unit vector with the same direction, or `None` for a
    /// zero (or non-finite) vector.
    #[inline]
    #[must_use]
    pub fn try_normalize(self) -> Option<Self> {
        let n = self.norm();
        (n > 0.0 && n.is_finite()).then(|| self / n)
    }

    /// Returns the unit vector with the same direction.
    ///
    /// A zero vector yields NaN components; use [`Vec3::try_normalize`] when
    /// the input may be degenerate.
    #[inline]
    #[must_use]
    pub fn normalize(self) -> Self {
        self / self.norm()
    }

    /// Euclidean distance to another point.
    #[inline]
    #[must_use]
    pub fn distance(self, o: Self) -> f64 {
        (self - o).norm()
    }

    /// Angle between two vectors in radians, in `[0, π]`.
    ///
    /// Uses `atan2(|a×b|, a·b)`, which stays accurate for nearly parallel
    /// and nearly anti-parallel vectors (unlike `acos`).
    #[inline]
    #[must_use]
    pub fn angle_to(self, o: Self) -> f64 {
        self.cross(o).norm().atan2(self.dot(o))
    }

    /// True when every component is finite.
    #[inline]
    #[must_use]
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    /// Any unit vector orthogonal to `self` (which need not be normalised).
    #[inline]
    #[must_use]
    pub fn any_orthogonal(self) -> Self {
        // Cross with the axis least aligned with `self` for good conditioning.
        let a = if self.x.abs() <= self.y.abs() && self.x.abs() <= self.z.abs() {
            Self::X
        } else if self.y.abs() <= self.z.abs() {
            Self::Y
        } else {
            Self::Z
        };
        self.cross(a).normalize()
    }

    /// Rotates `self` around the unit `axis` by `angle` radians (Rodrigues).
    #[inline]
    #[must_use]
    pub fn rotate_about(self, axis: Self, angle: f64) -> Self {
        let (s, c) = angle.sin_cos();
        self * c + axis.cross(self) * s + axis * (axis.dot(self) * (1.0 - c))
    }
}

impl Add for Vec3 {
    type Output = Self;
    #[inline]
    fn add(self, o: Self) -> Self {
        Self::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }
}

impl AddAssign for Vec3 {
    #[inline]
    fn add_assign(&mut self, o: Self) {
        *self = *self + o;
    }
}

impl Sub for Vec3 {
    type Output = Self;
    #[inline]
    fn sub(self, o: Self) -> Self {
        Self::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }
}

impl SubAssign for Vec3 {
    #[inline]
    fn sub_assign(&mut self, o: Self) {
        *self = *self - o;
    }
}

impl Mul<f64> for Vec3 {
    type Output = Self;
    #[inline]
    fn mul(self, s: f64) -> Self {
        Self::new(self.x * s, self.y * s, self.z * s)
    }
}

impl Mul<Vec3> for f64 {
    type Output = Vec3;
    #[inline]
    fn mul(self, v: Vec3) -> Vec3 {
        v * self
    }
}

impl Div<f64> for Vec3 {
    type Output = Self;
    #[inline]
    fn div(self, s: f64) -> Self {
        Self::new(self.x / s, self.y / s, self.z / s)
    }
}

impl Neg for Vec3 {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

/// A 3×3 matrix stored as rows. Rotation matrices built by this crate are
/// orthonormal, so the inverse is the transpose.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Mat3 {
    /// Rows of the matrix.
    pub rows: [Vec3; 3],
}

impl Mat3 {
    /// Identity matrix.
    pub const IDENTITY: Self = Self {
        rows: [Vec3::X, Vec3::Y, Vec3::Z],
    };

    /// Builds a matrix from its rows.
    #[inline]
    #[must_use]
    pub const fn from_rows(r0: Vec3, r1: Vec3, r2: Vec3) -> Self {
        Self { rows: [r0, r1, r2] }
    }

    /// Builds a matrix from its columns.
    #[inline]
    #[must_use]
    pub const fn from_cols(c0: Vec3, c1: Vec3, c2: Vec3) -> Self {
        Self::from_rows(
            Vec3::new(c0.x, c1.x, c2.x),
            Vec3::new(c0.y, c1.y, c2.y),
            Vec3::new(c0.z, c1.z, c2.z),
        )
    }

    /// Matrix–vector product.
    #[inline]
    #[must_use]
    pub fn mul_vec(&self, v: Vec3) -> Vec3 {
        Vec3::new(self.rows[0].dot(v), self.rows[1].dot(v), self.rows[2].dot(v))
    }

    /// Transposed-matrix–vector product (the inverse for a rotation).
    #[inline]
    #[must_use]
    pub fn mul_vec_t(&self, v: Vec3) -> Vec3 {
        self.rows[0] * v.x + self.rows[1] * v.y + self.rows[2] * v.z
    }

    /// Transpose.
    #[inline]
    #[must_use]
    pub fn transpose(&self) -> Self {
        Self::from_cols(self.rows[0], self.rows[1], self.rows[2])
    }

    /// Matrix product `self * o`.
    #[inline]
    #[must_use]
    pub fn mul_mat(&self, o: &Self) -> Self {
        let t = o.transpose();
        Self::from_rows(
            Vec3::new(
                self.rows[0].dot(t.rows[0]),
                self.rows[0].dot(t.rows[1]),
                self.rows[0].dot(t.rows[2]),
            ),
            Vec3::new(
                self.rows[1].dot(t.rows[0]),
                self.rows[1].dot(t.rows[1]),
                self.rows[1].dot(t.rows[2]),
            ),
            Vec3::new(
                self.rows[2].dot(t.rows[0]),
                self.rows[2].dot(t.rows[1]),
                self.rows[2].dot(t.rows[2]),
            ),
        )
    }

    /// Rotation matrix about the X axis (active rotation, right-handed).
    #[must_use]
    pub fn rot_x(a: f64) -> Self {
        let (s, c) = a.sin_cos();
        Self::from_rows(Vec3::X, Vec3::new(0.0, c, -s), Vec3::new(0.0, s, c))
    }

    /// Rotation matrix about the Y axis (active rotation, right-handed).
    #[must_use]
    pub fn rot_y(a: f64) -> Self {
        let (s, c) = a.sin_cos();
        Self::from_rows(Vec3::new(c, 0.0, s), Vec3::Y, Vec3::new(-s, 0.0, c))
    }

    /// Rotation matrix about the Z axis (active rotation, right-handed).
    #[must_use]
    pub fn rot_z(a: f64) -> Self {
        let (s, c) = a.sin_cos();
        Self::from_rows(Vec3::new(c, -s, 0.0), Vec3::new(s, c, 0.0), Vec3::Z)
    }
}

/// Bisection for the boundary of a predicate that is true at `lo` and false
/// somewhere towards `hi`: returns the last point found where it is true,
/// after `steps` halvings (so within `(hi − lo) / 2^steps` of the boundary).
///
/// The predicate is assumed monotone (true then false) on `[lo, hi]`.
pub fn bisect_last_true<P: FnMut(f64) -> bool>(mut lo: f64, mut hi: f64, steps: usize, mut pred: P) -> f64 {
    for _ in 0..steps {
        let mid = 0.5 * (lo + hi);
        if pred(mid) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo
}

/// Number of halvings needed to shrink an interval of `width` below `tol`.
#[must_use]
pub fn bisect_steps_for(width: f64, tol: f64) -> usize {
    if width <= tol {
        0
    } else {
        (width / tol).log2().ceil() as usize
    }
}

/// Wraps an angle in radians to `(-π, π]`.
#[inline]
#[must_use]
pub fn wrap_pi(a: f64) -> f64 {
    use std::f64::consts::{PI, TAU};
    let mut r = a % TAU;
    if r <= -PI {
        r += TAU;
    } else if r > PI {
        r -= TAU;
    }
    r
}

/// Wraps an angle in radians to `[0, 2π)`.
#[inline]
#[must_use]
pub fn wrap_two_pi(a: f64) -> f64 {
    use std::f64::consts::TAU;
    let r = a.rem_euclid(TAU);
    // rem_euclid can return TAU itself for tiny negative inputs.
    if r >= TAU { 0.0 } else { r }
}

/// Wraps a longitude in degrees to `[-180, 180)`; `180` maps to `-180`.
///
/// Use [`wrap_lon_deg_closed`] if `+180` must be preserved.
#[inline]
#[must_use]
pub fn wrap_lon_deg(lon: f64) -> f64 {
    let r = (lon + 180.0).rem_euclid(360.0) - 180.0;
    if r >= 180.0 { -180.0 } else { r }
}

/// Wraps a longitude difference in degrees to `(-180, 180]`: the shortest
/// signed step, with an exact half turn kept as `+180`.
#[inline]
#[must_use]
pub fn wrap_lon_delta_deg(d: f64) -> f64 {
    let r = (d + 180.0).rem_euclid(360.0) - 180.0;
    if r <= -180.0 { 180.0 } else { r }
}

/// Wraps a longitude in degrees to `(-180, 180]` (keeps `+180`).
#[inline]
#[must_use]
pub fn wrap_lon_deg_closed(lon: f64) -> f64 {
    let r = wrap_lon_deg(lon);
    if r == -180.0 && lon > 0.0 { 180.0 } else { r }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use std::f64::consts::{FRAC_PI_2, PI};

    #[test]
    fn cross_and_dot() {
        assert_eq!(Vec3::X.cross(Vec3::Y), Vec3::Z);
        assert_eq!(Vec3::X.dot(Vec3::Y), 0.0);
        assert_abs_diff_eq!(Vec3::new(3.0, 4.0, 0.0).norm(), 5.0);
    }

    #[test]
    fn normalize_zero_is_none() {
        assert!(Vec3::ZERO.try_normalize().is_none());
        assert!(Vec3::new(f64::INFINITY, 0.0, 0.0).try_normalize().is_none());
        assert_eq!(Vec3::new(0.0, 2.0, 0.0).try_normalize(), Some(Vec3::Y));
    }

    #[test]
    fn angle_to_is_stable() {
        assert_abs_diff_eq!(Vec3::X.angle_to(Vec3::Y), FRAC_PI_2);
        assert_abs_diff_eq!(Vec3::X.angle_to(-Vec3::X), PI);
        let tiny = Vec3::new(1.0, 1e-12, 0.0);
        assert_abs_diff_eq!(Vec3::X.angle_to(tiny), 1e-12, epsilon = 1e-20);
    }

    #[test]
    fn any_orthogonal_is_orthogonal() {
        for v in [Vec3::X, Vec3::new(1.0, 2.0, 3.0), Vec3::new(0.0, 0.0, -5.0)] {
            let o = v.any_orthogonal();
            assert_abs_diff_eq!(o.dot(v), 0.0, epsilon = 1e-12);
            assert_abs_diff_eq!(o.norm(), 1.0, epsilon = 1e-12);
        }
    }

    #[test]
    fn rotations() {
        let r = Mat3::rot_z(FRAC_PI_2);
        let v = r.mul_vec(Vec3::X);
        assert_abs_diff_eq!(v.y, 1.0, epsilon = 1e-15);
        let back = r.mul_vec_t(v);
        assert_abs_diff_eq!(back.x, 1.0, epsilon = 1e-15);
        let rr = Vec3::X.rotate_about(Vec3::Z, FRAC_PI_2);
        assert_abs_diff_eq!(rr.y, 1.0, epsilon = 1e-15);
        let m = Mat3::rot_x(0.3).mul_mat(&Mat3::rot_x(-0.3));
        assert_abs_diff_eq!(m.rows[1].y, 1.0, epsilon = 1e-15);
        assert_eq!(Mat3::IDENTITY.transpose(), Mat3::IDENTITY);
        let ry = Mat3::rot_y(FRAC_PI_2).mul_vec(Vec3::Z);
        assert_abs_diff_eq!(ry.x, 1.0, epsilon = 1e-15);
    }

    #[test]
    fn wrapping() {
        assert_abs_diff_eq!(wrap_pi(3.0 * PI), PI, epsilon = 1e-12);
        assert_abs_diff_eq!(wrap_pi(-PI), PI, epsilon = 1e-12);
        assert_abs_diff_eq!(wrap_two_pi(-0.5), 2.0 * PI - 0.5, epsilon = 1e-12);
        assert_eq!(wrap_two_pi(-1e-300), 0.0);
        assert_eq!(wrap_lon_deg(180.0), -180.0);
        assert_eq!(wrap_lon_deg(190.0), -170.0);
        assert_eq!(wrap_lon_deg(-190.0), 170.0);
        assert_eq!(wrap_lon_deg_closed(180.0), 180.0);
        assert_eq!(wrap_lon_deg_closed(-180.0), -180.0);
        assert_eq!(wrap_lon_deg_closed(540.0), 180.0);
    }

    #[test]
    fn sub_assign_and_scalar_times_vector() {
        let mut v = Vec3::new(3.0, 2.0, 1.0);
        v -= Vec3::new(1.0, 1.0, 1.0);
        assert_eq!(v, Vec3::new(2.0, 1.0, 0.0));
        assert_eq!(2.0 * Vec3::X, Vec3::new(2.0, 0.0, 0.0));
        let mut w = Vec3::ZERO;
        w += Vec3::Y;
        assert_eq!(w, Vec3::Y);
        assert_eq!(-Vec3::Z, Vec3::new(0.0, 0.0, -1.0));
    }

    #[test]
    fn wrap_lon_deg_maps_rounding_to_plus_180_onto_minus_180() {
        // One ulp below -180 (a smaller offset rounds back to exactly -180):
        // lon + 180 is a tiny negative number and rem_euclid(360) rounds it to
        // 360, which must still map to -180.
        let lon: f64 = -180.0 - 3e-14;
        assert!(lon < -180.0);
        assert_eq!((lon + 180.0).rem_euclid(360.0), 360.0);
        assert_eq!(wrap_lon_deg(lon), -180.0);
    }

    #[test]
    fn wrap_lon_delta_keeps_half_turn_positive_and_takes_shortest_step() {
        assert_eq!(wrap_lon_delta_deg(180.0), 180.0);
        assert_eq!(wrap_lon_delta_deg(-180.0), 180.0);
        assert_eq!(wrap_lon_delta_deg(190.0), -170.0);
        assert_eq!(wrap_lon_delta_deg(-190.0), 170.0);
        assert_eq!(wrap_lon_delta_deg(0.0), 0.0);
        assert_eq!(wrap_lon_delta_deg(540.0), 180.0);
    }

    #[test]
    fn bisect_finds_boundary_of_monotone_predicate() {
        let x = bisect_last_true(0.0, 10.0, 60, |t| t <= std::f64::consts::E);
        assert_abs_diff_eq!(x, std::f64::consts::E, epsilon = 1e-12);
        // Zero steps returns the lower bound unchanged.
        assert_eq!(bisect_last_true(1.0, 2.0, 0, |_| true), 1.0);
    }

    #[test]
    fn bisect_steps_for_width_and_tolerance() {
        assert_eq!(bisect_steps_for(1.0, 2.0), 0);
        assert_eq!(bisect_steps_for(1.0, 1.0), 0);
        assert_eq!(bisect_steps_for(1.0, 0.25), 2);
        assert_eq!(bisect_steps_for(60.0, 0.1), 10);
    }
}
