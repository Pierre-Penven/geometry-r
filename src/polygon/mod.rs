//! Spherical polygons with great-circle edges.
//!
//! # Model
//! * Vertices are unit vectors built from `(lon, lat)` in degrees
//!   ([`LonLat`], GeoJSON order). Geodetic latitude is treated as spherical,
//!   as GeoJSON and S2 do.
//! * Edges are **great-circle minor arcs**. There is no special meaning to
//!   ±180° longitude: the antimeridian and the poles are ordinary places.
//! * Every ring bounds the region **on its left** seen from outside the
//!   sphere. Exterior rings are counter-clockwise, holes clockwise. The
//!   polygon is the intersection of the left regions of its rings.
//!
//! Boolean operations live in [`boolean`]; RFC 7946 output and bounding
//! boxes in [`antimeridian`].

pub mod antimeridian;
pub mod boolean;
mod planar;

use crate::ellipsoid::{EARTH_AUTHALIC_RADIUS_M, Ellipsoid};
use crate::error::{GeoError, Result};
use crate::geodesic::slerp;
use crate::math::Vec3;
use std::f64::consts::{FRAC_PI_2, PI, TAU};

/// A `[longitude, latitude]` pair in degrees (GeoJSON order).
pub type LonLat = [f64; 2];

/// Two vertices closer than this (chord on the unit sphere, ≈ 6 µm on the
/// Earth) are merged.
const DUP_EPS2: f64 = 1e-24;

/// `[lon, lat]` degrees → unit vector.
#[inline]
#[must_use]
pub fn lonlat_to_unit(ll: LonLat) -> Vec3 {
    Ellipsoid::normal(ll[1].to_radians(), ll[0].to_radians())
}

/// Unit vector → `[lon, lat]` degrees, longitude in `(-180, 180]`.
#[inline]
#[must_use]
pub fn unit_to_lonlat(v: Vec3) -> LonLat {
    [v.y.atan2(v.x).to_degrees(), v.z.atan2(v.x.hypot(v.y)).to_degrees()]
}

/// A spherical cap: all points within `radius_rad` of `center`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cap {
    center: Vec3,
    radius_rad: f64,
    cos_radius: f64,
}

impl Cap {
    /// Builds a cap. The center is normalised; the radius is clamped to
    /// `[0, π]`.
    ///
    /// # Errors
    /// [`GeoError::InvalidArgument`] for a zero center or a non-finite radius.
    pub fn new(center: Vec3, radius_rad: f64) -> Result<Self> {
        let c = center
            .try_normalize()
            .ok_or_else(|| GeoError::InvalidArgument("zero cap center".into()))?;
        if !radius_rad.is_finite() {
            return Err(GeoError::InvalidArgument("non-finite cap radius".into()));
        }
        Ok(Self::new_unchecked(c, radius_rad))
    }

    fn new_unchecked(center: Vec3, radius_rad: f64) -> Self {
        let r = radius_rad.clamp(0.0, PI);
        Self {
            center,
            radius_rad: r,
            cos_radius: if r >= PI { -1.0 } else { r.cos() },
        }
    }

    /// The whole sphere.
    #[must_use]
    pub fn full() -> Self {
        Self::new_unchecked(Vec3::Z, PI)
    }

    /// Unit center.
    #[must_use]
    pub fn center(&self) -> Vec3 {
        self.center
    }

    /// Angular radius in radians.
    #[must_use]
    pub fn radius_rad(&self) -> f64 {
        self.radius_rad
    }

    /// True for the whole sphere.
    #[must_use]
    pub fn is_full(&self) -> bool {
        self.radius_rad >= PI
    }

    /// True when the unit vector `p` is inside the cap (with a tiny
    /// tolerance so that boundary vertices are included).
    #[inline]
    #[must_use]
    pub fn contains(&self, p: Vec3) -> bool {
        self.is_full() || p.dot(self.center) >= self.cos_radius - 1e-12
    }

    /// True when the two caps overlap.
    #[must_use]
    pub fn intersects(&self, o: &Cap) -> bool {
        self.is_full() || o.is_full() || self.center.angle_to(o.center) <= self.radius_rad + o.radius_rad + 1e-12
    }

    /// Smallest cap containing both caps.
    #[must_use]
    pub fn union(&self, o: &Cap) -> Cap {
        if self.is_full() || o.is_full() {
            return Cap::full();
        }
        let d = self.center.angle_to(o.center);
        if d + o.radius_rad <= self.radius_rad {
            return *self;
        }
        if d + self.radius_rad <= o.radius_rad {
            return *o;
        }
        let r = 0.5 * (d + self.radius_rad + o.radius_rad);
        if r >= PI {
            return Cap::full();
        }
        // Move from self.center toward o.center by (r - self.radius).
        let axis = self.center.cross(o.center);
        let center = match axis.try_normalize() {
            Some(ax) => self.center.rotate_about(ax, r - self.radius_rad),
            None => return Cap::full(), // antipodal centers
        };
        Cap::new_unchecked(center, r + 1e-12)
    }
}

/// How ring orientation is decided when building a polygon.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Orientation {
    /// Ignore the input winding: the exterior bounds the smaller of the two
    /// regions it defines, each hole cuts out the smaller region it defines.
    /// Safe default for user input (cannot represent a polygon larger than a
    /// hemisphere).
    #[default]
    SmallestArea,
    /// Trust the input: exterior counter-clockwise, holes clockwise (seen
    /// from outside the sphere). Needed for regions larger than a hemisphere.
    AsGiven,
}

/// A spherical polygon: one exterior ring and zero or more holes.
#[derive(Debug, Clone, PartialEq)]
pub struct Polygon {
    rings: Vec<Vec<Vec3>>,
    cap: Cap,
    /// True when the polygon fits in a cap smaller than a hemisphere.
    small: bool,
    /// A point known to be inside the polygon.
    inside_ref: Vec3,
}

impl Polygon {
    /// Builds a polygon from lon/lat rings (first ring exterior, the rest
    /// holes). Rings may be closed or open; duplicates and spikes are
    /// removed. Uses [`Orientation::SmallestArea`].
    ///
    /// # Errors
    /// [`GeoError::InvalidCoordinate`] for non-finite or out-of-range
    /// coordinates; [`GeoError::DegeneratePolygon`] for a ring with fewer
    /// than three distinct vertices or no rings.
    pub fn from_lonlat(rings: &[Vec<LonLat>]) -> Result<Self> {
        Self::from_lonlat_with(rings, Orientation::SmallestArea)
    }

    /// Like [`Polygon::from_lonlat`] with an explicit orientation policy.
    ///
    /// # Errors
    /// See [`Polygon::from_lonlat`].
    pub fn from_lonlat_with(rings: &[Vec<LonLat>], orientation: Orientation) -> Result<Self> {
        let mut out = Vec::with_capacity(rings.len());
        for ring in rings {
            let mut r = Vec::with_capacity(ring.len());
            for &[lon, lat] in ring {
                if !(lon.is_finite() && lat.is_finite()) || lat.abs() > 90.0 {
                    return Err(GeoError::InvalidCoordinate(format!("[{lon}, {lat}]")));
                }
                r.push(lonlat_to_unit([lon, lat]));
            }
            out.push(r);
        }
        Self::from_unit_rings(out, orientation)
    }

    /// Builds a polygon from rings of unit vectors.
    ///
    /// # Errors
    /// [`GeoError::DegeneratePolygon`] for no rings or a ring with fewer
    /// than three distinct vertices; [`GeoError::InvalidCoordinate`] for
    /// non-finite vectors.
    pub fn from_unit_rings(rings: Vec<Vec<Vec3>>, orientation: Orientation) -> Result<Self> {
        if rings.is_empty() {
            return Err(GeoError::DegeneratePolygon("no rings".into()));
        }
        let mut cleaned = Vec::with_capacity(rings.len());
        for (i, ring) in rings.into_iter().enumerate() {
            let mut r = Vec::with_capacity(ring.len());
            for v in ring {
                let v = v
                    .try_normalize()
                    .ok_or_else(|| GeoError::InvalidCoordinate("zero or non-finite vertex".into()))?;
                r.push(v);
            }
            let mut r = clean_ring(r);
            if r.len() < 3 {
                return Err(GeoError::DegeneratePolygon(format!(
                    "ring {i} has fewer than 3 distinct vertices"
                )));
            }
            if orientation == Orientation::SmallestArea {
                let s = ring_signed_area(&r);
                if (i == 0 && s < 0.0) || (i > 0 && s > 0.0) {
                    r.reverse();
                }
            }
            cleaned.push(r);
        }
        Ok(Self::build(cleaned))
    }

    /// Internal constructor for already-clean, correctly oriented rings.
    pub(crate) fn build(rings: Vec<Vec<Vec3>>) -> Self {
        let (cap, small) = bounding_cap(&rings[0]);
        let inside_ref = inside_reference(&rings);
        Self {
            rings,
            cap,
            small,
            inside_ref,
        }
    }

    /// A circle (cap boundary) of `n ≥ 3` vertices around `center` with
    /// angular radius `radius_rad ∈ (0, π)`. Radii above 90° produce
    /// polygons larger than a hemisphere.
    ///
    /// # Errors
    /// [`GeoError::InvalidArgument`] for `n < 3` or a radius outside `(0, π)`.
    pub fn circle(center: Vec3, radius_rad: f64, n: usize) -> Result<Self> {
        if n < 3 {
            return Err(GeoError::InvalidArgument("circle needs at least 3 vertices".into()));
        }
        if !(radius_rad > 0.0 && radius_rad < PI) {
            return Err(GeoError::InvalidArgument("circle radius must be in (0, π)".into()));
        }
        let c = center
            .try_normalize()
            .ok_or_else(|| GeoError::InvalidArgument("zero circle center".into()))?;
        let start = c.rotate_about(c.any_orthogonal(), radius_rad);
        // Rotating counter-clockwise about the outward center keeps the cap on the left.
        let ring = (0..n)
            .map(|i| start.rotate_about(c, TAU * i as f64 / n as f64))
            .collect();
        Self::from_unit_rings(vec![ring], Orientation::AsGiven)
    }

    /// All rings; index 0 is the exterior.
    #[must_use]
    pub fn rings(&self) -> &[Vec<Vec3>] {
        &self.rings
    }

    /// Exterior ring (counter-clockwise, open).
    #[must_use]
    pub fn exterior(&self) -> &[Vec3] {
        &self.rings[0]
    }

    /// Hole rings (clockwise, open).
    #[must_use]
    pub fn holes(&self) -> &[Vec<Vec3>] {
        &self.rings[1..]
    }

    /// Bounding cap.
    #[must_use]
    pub fn cap(&self) -> Cap {
        self.cap
    }

    /// True when the polygon lies within a cap smaller than a hemisphere.
    #[must_use]
    pub fn is_small(&self) -> bool {
        self.small
    }

    /// Area in steradians (unit sphere).
    #[must_use]
    pub fn area_sr(&self) -> f64 {
        let ext = ring_signed_area(&self.rings[0]).rem_euclid(4.0 * PI);
        let holes: f64 = self.rings[1..]
            .iter()
            .map(|h| (-ring_signed_area(h)).rem_euclid(4.0 * PI))
            .sum();
        (ext - holes).max(0.0)
    }

    /// Area in m² on a sphere of radius `radius_m`.
    #[must_use]
    pub fn area_m2_with_radius(&self, radius_m: f64) -> f64 {
        self.area_sr() * radius_m * radius_m
    }

    /// Area in m² on the authalic Earth sphere (close to the WGS-84
    /// ellipsoidal area for large regions).
    #[must_use]
    pub fn area_m2(&self) -> f64 {
        self.area_m2_with_radius(EARTH_AUTHALIC_RADIUS_M)
    }

    /// Perimeter (all rings) in radians on the unit sphere.
    #[must_use]
    pub fn perimeter_rad(&self) -> f64 {
        self.rings
            .iter()
            .map(|r| ring_edges(r).map(|(a, b)| a.angle_to(b)).sum::<f64>())
            .sum()
    }

    /// Normalised vector sum of the area-weighted triangle centroids: a good
    /// "center" for small polygons. Returns the cap center for degenerate
    /// cases.
    #[must_use]
    pub fn centroid(&self) -> Vec3 {
        let mut acc = Vec3::ZERO;
        // Hole rings are clockwise, so their triangles subtract themselves.
        for ring in &self.rings {
            let c = ring_center(ring);
            for (a, b) in ring_edges(ring) {
                let area = triangle_signed_area(c, a, b);
                acc += (c + a + b).normalize() * area;
            }
        }
        acc.try_normalize().unwrap_or(self.cap.center)
    }

    /// True when the unit vector `p` is inside the polygon (boundary points
    /// may go either way).
    #[must_use]
    pub fn contains(&self, p: Vec3) -> bool {
        if !self.cap.contains(p) {
            return false;
        }
        let r = self.inside_ref;
        // Path p → r; split when nearly antipodal so both legs are minor arcs.
        let crossings = if p.dot(r) < -0.9 {
            let m = slerp(p, r, 0.5);
            self.count_crossings(p, m) + self.count_crossings(m, r)
        } else {
            self.count_crossings(p, r)
        };
        crossings % 2 == 0
    }

    /// [`Polygon::contains`] for a `[lon, lat]` point in degrees.
    #[must_use]
    pub fn contains_lonlat(&self, ll: LonLat) -> bool {
        self.contains(lonlat_to_unit(ll))
    }

    fn count_crossings(&self, a: Vec3, b: Vec3) -> usize {
        self.rings
            .iter()
            .map(|ring| ring_edges(ring).filter(|&(c, d)| arcs_cross(a, b, c, d)).count())
            .sum()
    }

    /// Inserts vertices so that no edge is longer than `max_edge_rad`.
    ///
    /// # Errors
    /// [`GeoError::InvalidArgument`] for a non-positive limit.
    pub fn densify(&self, max_edge_rad: f64) -> Result<Self> {
        if max_edge_rad.is_nan() || max_edge_rad <= 0.0 {
            return Err(GeoError::InvalidArgument("max edge must be > 0".into()));
        }
        let rings = self
            .rings
            .iter()
            .map(|ring| {
                let mut out = Vec::with_capacity(ring.len());
                for (a, b) in ring_edges(ring) {
                    let n = (a.angle_to(b) / max_edge_rad).ceil().max(1.0) as usize;
                    for i in 0..n {
                        out.push(slerp(a, b, i as f64 / n as f64));
                    }
                }
                out
            })
            .collect();
        Ok(Self::build(rings))
    }

    /// Rings as closed `[lon, lat]` sequences (first vertex repeated),
    /// without any antimeridian splitting. See
    /// [`antimeridian::to_rfc7946`] for standard-compliant output.
    #[must_use]
    pub fn to_lonlat_rings(&self) -> Vec<Vec<LonLat>> {
        self.rings
            .iter()
            .map(|r| {
                let mut v: Vec<LonLat> = r.iter().map(|&p| unit_to_lonlat(p)).collect();
                v.push(v[0]);
                v
            })
            .collect()
    }
}

/// A set of non-overlapping polygons.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct MultiPolygon {
    /// Member polygons.
    pub polygons: Vec<Polygon>,
}

impl From<Polygon> for MultiPolygon {
    fn from(p: Polygon) -> Self {
        Self { polygons: vec![p] }
    }
}

impl MultiPolygon {
    /// Builds from polygons (assumed non-overlapping; call
    /// [`MultiPolygon::dissolve`] otherwise).
    #[must_use]
    pub fn new(polygons: Vec<Polygon>) -> Self {
        Self { polygons }
    }

    /// True when there is no polygon.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.polygons.is_empty()
    }

    /// Total area in steradians.
    #[must_use]
    pub fn area_sr(&self) -> f64 {
        self.polygons.iter().map(Polygon::area_sr).sum()
    }

    /// Total area in m² on the authalic Earth sphere.
    #[must_use]
    pub fn area_m2(&self) -> f64 {
        self.polygons.iter().map(Polygon::area_m2).sum()
    }

    /// True when any member contains `p`.
    #[must_use]
    pub fn contains(&self, p: Vec3) -> bool {
        self.polygons.iter().any(|q| q.contains(p))
    }

    /// [`MultiPolygon::contains`] for a `[lon, lat]` point in degrees.
    #[must_use]
    pub fn contains_lonlat(&self, ll: LonLat) -> bool {
        self.contains(lonlat_to_unit(ll))
    }

    /// Bounding cap of all members (`None` when empty).
    #[must_use]
    pub fn cap(&self) -> Option<Cap> {
        let mut it = self.polygons.iter();
        let first = it.next()?.cap();
        Some(it.fold(first, |acc, p| acc.union(&p.cap())))
    }
}

/// Iterator over the edges `(v[i], v[i+1])` of a closed ring.
#[inline]
pub(crate) fn ring_edges(ring: &[Vec3]) -> impl Iterator<Item = (Vec3, Vec3)> + '_ {
    let n = ring.len();
    (0..n).map(move |i| (ring[i], ring[(i + 1) % n]))
}

/// Signed area of the spherical triangle `abc` (positive counter-clockwise).
#[inline]
pub(crate) fn triangle_signed_area(a: Vec3, b: Vec3, c: Vec3) -> f64 {
    let num = a.dot(b.cross(c));
    let den = 1.0 + a.dot(b) + b.dot(c) + c.dot(a);
    2.0 * num.atan2(den)
}

fn ring_center(ring: &[Vec3]) -> Vec3 {
    let s = ring.iter().fold(Vec3::ZERO, |acc, &v| acc + v);
    s.try_normalize().filter(|_| s.norm() > 1e-9).unwrap_or(ring[0])
}

/// Signed area (steradians) of the region left of the ring, as a fan from
/// the ring center. Positive and below 2π for a small counter-clockwise
/// ring; negative for a small clockwise ring. For a ring whose left side is
/// larger than a hemisphere the true area is `value.rem_euclid(4π)`.
#[must_use]
pub fn ring_signed_area(ring: &[Vec3]) -> f64 {
    let c = ring_center(ring);
    ring_edges(ring).map(|(a, b)| triangle_signed_area(c, a, b)).sum()
}

/// True when minor arcs `ab` and `cd` cross.
///
/// A vertex lying exactly on the great circle of `ab` is treated as lying on
/// its negative side (half-open rule), so a path passing exactly through a
/// ring vertex is counted once, keeping crossing parity correct.
#[inline]
pub(crate) fn arcs_cross(a: Vec3, b: Vec3, c: Vec3, d: Vec3) -> bool {
    let n1 = a.cross(b);
    if (n1.dot(c) > 0.0) == (n1.dot(d) > 0.0) {
        return false;
    }
    let n2 = c.cross(d);
    if (n2.dot(a) > 0.0) == (n2.dot(b) > 0.0) {
        return false;
    }
    // The circles meet at ±x; keep the point on the minor arc ab and check
    // that it is also on the minor arc cd (not its antipode).
    let mut x = n1.cross(n2);
    if x.dot(a + b) < 0.0 {
        x = -x;
    }
    x.dot(c + d) > 0.0
}

/// Removes closing duplicates, consecutive duplicates and spikes.
fn clean_ring(ring: Vec<Vec3>) -> Vec<Vec3> {
    let same = |a: Vec3, b: Vec3| (a - b).norm2() < DUP_EPS2;
    let mut out: Vec<Vec3> = Vec::with_capacity(ring.len());
    for p in ring {
        if out.last().is_some_and(|&l| same(l, p)) {
            continue;
        }
        if out.len() >= 2 && same(out[out.len() - 2], p) {
            out.pop(); // a b a → a
            continue;
        }
        out.push(p);
    }
    // Cyclic clean-up (closing point, spikes across the seam).
    loop {
        let n = out.len();
        if n < 3 {
            break;
        }
        if same(out[0], out[n - 1]) {
            out.pop();
            continue;
        }
        if same(out[n - 2], out[0]) {
            out.truncate(n - 1); // ... a, b | a ... : drop spike tip b
            out.pop();
            continue;
        }
        if same(out[n - 1], out[1]) {
            out.remove(0); // spike tip at index 0
            out.pop();
            continue;
        }
        break;
    }
    out
}

/// Bounding cap of the region left of the exterior ring, and whether the
/// polygon is "small" (inside a cap of radius < 90°).
fn bounding_cap(exterior: &[Vec3]) -> (Cap, bool) {
    let s = ring_signed_area(exterior);
    if !(s > 0.0 && s <= TAU) {
        return (Cap::full(), false);
    }
    let (c, r) = enclosing_cap(exterior.iter().copied(), ring_center(exterior));
    if r < FRAC_PI_2 - 1e-9 {
        (Cap::new_unchecked(c, r + 1e-12), true)
    } else {
        (Cap::full(), false)
    }
}

/// Near-minimal cap enclosing a point set: Bădoiu–Clarkson iterations
/// (repeatedly step towards the farthest point) from `start`, keeping the
/// best center found. Returns `(center, radius)`.
pub(crate) fn enclosing_cap<I>(points: I, start: Vec3) -> (Vec3, f64)
where
    I: Iterator<Item = Vec3> + Clone,
{
    let farthest = |c: Vec3| {
        points
            .clone()
            .map(|v| (v, c.angle_to(v)))
            .fold((c, 0.0), |acc, x| if x.1 > acc.1 { x } else { acc })
    };
    let mut c = start;
    let (mut far, mut r) = farthest(c);
    let (mut best_c, mut best_r) = (c, r);
    for k in 0..48 {
        let step = 1.0 / f64::from(k + 2);
        match slerp_toward(c, far, r * step) {
            Some(nc) => c = nc,
            None => break,
        }
        (far, r) = farthest(c);
        if r < best_r {
            best_c = c;
            best_r = r;
        }
    }
    (best_c, best_r)
}

/// Rotates `a` towards `b` by `angle` radians (None when they coincide or
/// are antipodal).
fn slerp_toward(a: Vec3, b: Vec3, angle: f64) -> Option<Vec3> {
    let axis = a.cross(b).try_normalize()?;
    Some(a.rotate_about(axis, angle))
}

/// A point just left of the longest exterior edge: inside the polygon.
fn inside_reference(rings: &[Vec<Vec3>]) -> Vec3 {
    let ext = &rings[0];
    let (a, b) = ring_edges(ext)
        .max_by(|x, y| x.0.angle_to(x.1).total_cmp(&y.0.angle_to(y.1)))
        .unwrap_or((ext[0], ext[1]));
    let len = a.angle_to(b);
    let mid = slerp(a, b, 0.5);
    let left = a.cross(b).normalize();
    // Offset small relative to the edge, and below ~0.6 m on the Earth.
    let eps = (len * 1e-3).min(1e-7);
    (mid + left * eps).normalize()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    fn square(lon0: f64, lat0: f64, size: f64) -> Polygon {
        Polygon::from_lonlat(&[vec![
            [lon0, lat0],
            [lon0 + size, lat0],
            [lon0 + size, lat0 + size],
            [lon0, lat0 + size],
            [lon0, lat0],
        ]])
        .unwrap()
    }

    #[test]
    fn lonlat_round_trip() {
        let ll = [-179.5, 45.0];
        let back = unit_to_lonlat(lonlat_to_unit(ll));
        assert_abs_diff_eq!(back[0], ll[0], epsilon = 1e-12);
        assert_abs_diff_eq!(back[1], ll[1], epsilon = 1e-12);
    }

    #[test]
    fn construction_errors() {
        assert!(matches!(Polygon::from_lonlat(&[]), Err(GeoError::DegeneratePolygon(_))));
        assert!(matches!(
            Polygon::from_lonlat(&[vec![[0.0, 0.0], [1.0, 0.0], [0.0, 0.0]]]),
            Err(GeoError::DegeneratePolygon(_))
        ));
        assert!(matches!(
            Polygon::from_lonlat(&[vec![[0.0, 0.0], [1.0, 95.0], [0.0, 1.0]]]),
            Err(GeoError::InvalidCoordinate(_))
        ));
        assert!(matches!(
            Polygon::from_lonlat(&[vec![[f64::NAN, 0.0], [1.0, 0.0], [0.0, 1.0]]]),
            Err(GeoError::InvalidCoordinate(_))
        ));
        assert!(matches!(
            Polygon::from_unit_rings(vec![vec![Vec3::ZERO, Vec3::X, Vec3::Y]], Orientation::AsGiven),
            Err(GeoError::InvalidCoordinate(_))
        ));
        assert!(matches!(
            Polygon::circle(Vec3::Z, 0.1, 2),
            Err(GeoError::InvalidArgument(_))
        ));
        assert!(matches!(
            Polygon::circle(Vec3::Z, 4.0, 10),
            Err(GeoError::InvalidArgument(_))
        ));
        assert!(matches!(
            Polygon::circle(Vec3::ZERO, 0.1, 10),
            Err(GeoError::InvalidArgument(_))
        ));
    }

    #[test]
    fn ring_cleaning() {
        // Closing duplicate, consecutive duplicate and spike.
        let p = Polygon::from_lonlat(&[vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [2.5, 0.5],
            [2.0, 0.0],
            [2.0, 1.0],
            [0.0, 1.0],
            [0.0, 0.0],
        ]])
        .unwrap();
        assert_eq!(p.exterior().len(), 5);
        // Spike across the seam.
        let r = clean_ring(vec![Vec3::X, Vec3::Y, Vec3::Z, Vec3::Y]);
        assert_eq!(r.len(), 2);
    }

    #[test]
    fn orientation_is_normalised() {
        let ccw = square(10.0, 10.0, 5.0);
        let cw = Polygon::from_lonlat(&[vec![[10.0, 10.0], [10.0, 15.0], [15.0, 15.0], [15.0, 10.0]]]).unwrap();
        assert_abs_diff_eq!(ccw.area_sr(), cw.area_sr(), epsilon = 1e-15);
        assert!(ring_signed_area(cw.exterior()) > 0.0);
        assert!(ccw.is_small());
        // AsGiven keeps the big side.
        let big = Polygon::from_lonlat_with(
            &[vec![[10.0, 10.0], [10.0, 15.0], [15.0, 15.0], [15.0, 10.0]]],
            Orientation::AsGiven,
        )
        .unwrap();
        assert_abs_diff_eq!(big.area_sr(), 4.0 * PI - ccw.area_sr(), epsilon = 1e-12);
        assert!(!big.is_small());
        assert!(big.cap().is_full());
    }

    #[test]
    fn cap_area_matches_formula() {
        for &r in &[0.01, 0.5, 1.4, 2.0, 3.0] {
            let p = Polygon::circle(lonlat_to_unit([179.0, 80.0]), r, 2000).unwrap();
            let exact = TAU * (1.0 - r.cos());
            assert_abs_diff_eq!(p.area_sr() / exact, 1.0, epsilon = 2e-5);
            assert_eq!(p.is_small(), r < FRAC_PI_2);
        }
    }

    #[test]
    fn area_with_hole_and_earth_units() {
        let outer = vec![[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]];
        let hole = vec![[2.0, 2.0], [4.0, 2.0], [4.0, 4.0], [2.0, 4.0]];
        let p = Polygon::from_lonlat(&[outer.clone(), hole.clone()]).unwrap();
        let o = Polygon::from_lonlat(&[outer]).unwrap();
        let h = Polygon::from_lonlat(&[hole]).unwrap();
        assert_abs_diff_eq!(p.area_sr(), o.area_sr() - h.area_sr(), epsilon = 1e-14);
        // 1°×1° at the equator ≈ 12 364 km².
        let s = square(0.0, 0.0, 1.0);
        assert_abs_diff_eq!(s.area_m2() / 1e6, 12_364.0, epsilon = 15.0);
        assert_eq!(p.holes().len(), 1);
        assert!(p.perimeter_rad() > o.perimeter_rad());
    }

    #[test]
    fn contains_across_antimeridian_and_pole() {
        let p = Polygon::from_lonlat(&[vec![[170.0, -10.0], [-170.0, -10.0], [-170.0, 10.0], [170.0, 10.0]]]).unwrap();
        assert!(p.contains_lonlat([180.0, 0.0]));
        assert!(p.contains_lonlat([-175.0, 5.0]));
        assert!(!p.contains_lonlat([0.0, 0.0]));
        assert!(!p.contains_lonlat([160.0, 0.0]));
        assert_abs_diff_eq!(p.area_sr(), square(0.0, -10.0, 20.0).area_sr(), epsilon = 1e-12);
        // Polar cap ring.
        let ring: Vec<LonLat> = (0..36).map(|i| [-180.0 + 10.0 * f64::from(i), 70.0]).collect();
        let cap = Polygon::from_lonlat(&[ring]).unwrap();
        assert!(cap.contains_lonlat([0.0, 90.0]));
        assert!(cap.contains_lonlat([123.0, 75.0]));
        assert!(!cap.contains_lonlat([0.0, 60.0]));
        assert!(!cap.contains_lonlat([0.0, -90.0]));
    }

    #[test]
    fn contains_large_polygons() {
        // Circle of radius 120°: larger than a hemisphere.
        let c = lonlat_to_unit([0.0, 0.0]);
        let p = Polygon::circle(c, 120f64.to_radians(), 720).unwrap();
        assert!(p.contains(c));
        assert!(p.contains(lonlat_to_unit([100.0, 0.0])));
        assert!(p.contains(lonlat_to_unit([0.0, 90.0])));
        assert!(!p.contains(lonlat_to_unit([180.0, 0.0])));
        assert!(!p.contains(lonlat_to_unit([150.0, 10.0])));
        // Complement of a small square (AsGiven, clockwise exterior).
        let comp = Polygon::from_lonlat_with(
            &[vec![[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0]]],
            Orientation::AsGiven,
        )
        .unwrap();
        assert!(!comp.contains_lonlat([0.5, 0.5]));
        assert!(comp.contains_lonlat([-179.5, -0.5]));
        assert!(comp.contains_lonlat([0.5, -0.5]));
    }

    #[test]
    fn contains_with_hole() {
        let p = Polygon::from_lonlat(&[
            vec![[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]],
            vec![[2.0, 2.0], [4.0, 2.0], [4.0, 4.0], [2.0, 4.0]],
        ])
        .unwrap();
        assert!(p.contains_lonlat([1.0, 1.0]));
        assert!(!p.contains_lonlat([3.0, 3.0]));
        assert!(!p.contains_lonlat([11.0, 3.0]));
    }

    #[test]
    fn caps() {
        let a = Cap::new(Vec3::X, 0.1).unwrap();
        let b = Cap::new(Vec3::Y, 0.1).unwrap();
        assert!(!a.intersects(&b));
        let u = a.union(&b);
        assert!(u.contains(Vec3::X) && u.contains(Vec3::Y));
        assert_abs_diff_eq!(u.radius_rad(), FRAC_PI_4 + 0.1, epsilon = 1e-9);
        assert_eq!(a.union(&Cap::new(Vec3::X, 0.05).unwrap()), a);
        assert!(
            Cap::new(Vec3::X, 2.0)
                .unwrap()
                .union(&Cap::new(-Vec3::X, 2.0).unwrap())
                .is_full()
        );
        assert!(matches!(Cap::new(Vec3::ZERO, 1.0), Err(GeoError::InvalidArgument(_))));
        assert!(matches!(Cap::new(Vec3::X, f64::NAN), Err(GeoError::InvalidArgument(_))));
        assert!(Cap::full().contains(-Vec3::Z));
        let mp = MultiPolygon::new(vec![square(0.0, 0.0, 1.0), square(50.0, 0.0, 1.0)]);
        assert!(mp.cap().unwrap().contains(lonlat_to_unit([25.0, 0.5])));
        assert!(MultiPolygon::default().cap().is_none());
    }

    use std::f64::consts::FRAC_PI_4;

    #[test]
    fn densify_and_centroid() {
        let s = square(0.0, 0.0, 10.0);
        let d = s.densify(1f64.to_radians()).unwrap();
        assert!(d.exterior().len() >= 40);
        assert_abs_diff_eq!(d.area_sr(), s.area_sr(), epsilon = 1e-12);
        let c = unit_to_lonlat(s.centroid());
        assert_abs_diff_eq!(c[0], 5.0, epsilon = 0.01);
        assert_abs_diff_eq!(c[1], 5.0, epsilon = 0.1);
        assert!(matches!(s.densify(0.0), Err(GeoError::InvalidArgument(_))));
        let rings = s.to_lonlat_rings();
        assert_eq!(rings[0].first(), rings[0].last());
    }

    #[test]
    fn multipolygon_basics() {
        let mp: MultiPolygon = square(0.0, 0.0, 1.0).into();
        assert!(!mp.is_empty());
        assert!(mp.contains_lonlat([0.5, 0.5]));
        assert_abs_diff_eq!(mp.area_m2(), square(0.0, 0.0, 1.0).area_m2());
        assert_abs_diff_eq!(mp.area_sr(), square(0.0, 0.0, 1.0).area_sr());
    }

    #[test]
    fn arc_crossing() {
        let a = lonlat_to_unit([0.0, -1.0]);
        let b = lonlat_to_unit([0.0, 1.0]);
        let c = lonlat_to_unit([-1.0, 0.0]);
        let d = lonlat_to_unit([1.0, 0.0]);
        assert!(arcs_cross(a, b, c, d));
        assert!(!arcs_cross(a, b, -c, -d)); // antipodal arc: no crossing
        assert!(!arcs_cross(a, c, b, d));
    }

    #[test]
    fn cap_union_with_full_cap_is_full() {
        let a = Cap::new(Vec3::X, 0.1).unwrap();
        assert!(a.union(&Cap::full()).is_full());
        assert!(Cap::full().union(&a).is_full());
        assert!(Cap::full().intersects(&a));
    }

    #[test]
    fn cap_union_returns_containing_cap() {
        let small = Cap::new(Vec3::X, 0.05).unwrap();
        let big = Cap::new(Vec3::X, 0.5).unwrap();
        assert_eq!(small.union(&big), big);
        assert_eq!(big.union(&small), big);
    }

    #[test]
    fn cap_union_of_antipodal_caps_is_full() {
        let a = Cap::new(Vec3::X, 0.1).unwrap();
        let b = Cap::new(-Vec3::X, 0.1).unwrap();
        assert!(a.union(&b).is_full());
    }

    #[test]
    fn cap_accessors() {
        let c = Cap::new(Vec3::new(0.0, 0.0, 2.0), 0.3).unwrap();
        assert_eq!(c.center(), Vec3::Z);
        assert_eq!(c.radius_rad(), 0.3);
        assert!(!c.is_full());
        assert!(Cap::new(Vec3::Z, 10.0).unwrap().is_full(), "radius clamps to π");
    }

    #[test]
    fn clean_ring_removes_spike_closing_on_first_vertex() {
        let a = lonlat_to_unit([0.0, 0.0]);
        let b = lonlat_to_unit([1.0, 0.0]);
        let c = lonlat_to_unit([1.0, 1.0]);
        let d = lonlat_to_unit([-1.0, -1.0]);
        // A → D → A at the seam is a spike: D goes away.
        let r = clean_ring(vec![a, b, c, a, d]);
        assert_eq!(r, vec![a, b, c]);
    }

    #[test]
    fn long_thin_polygon_is_not_small_but_still_correct() {
        let p = Polygon::from_lonlat(&[vec![
            [0.0, 0.0],
            [100.0, 0.0],
            [200.0, 0.0],
            [200.0, 1.0],
            [100.0, 1.0],
            [0.0, 1.0],
        ]])
        .unwrap();
        assert!(!p.is_small(), "vertices span more than a hemisphere");
        assert!(p.cap().is_full());
        assert!(p.contains_lonlat([150.0, 0.5]));
        assert!(!p.contains_lonlat([150.0, 2.0]));
        assert!(!p.contains_lonlat([-90.0, 0.5]));
    }

    #[test]
    fn enclosing_cap_of_antipodal_points_stops_cleanly() {
        let (c, r) = enclosing_cap([Vec3::X, -Vec3::X].into_iter(), Vec3::X);
        assert_eq!(c, Vec3::X);
        assert_abs_diff_eq!(r, PI, epsilon = 1e-12);
    }

    #[test]
    fn enclosing_cap_converges_to_minimal_radius() {
        let pts = [
            lonlat_to_unit([-10.0, 0.0]),
            lonlat_to_unit([10.0, 0.0]),
            lonlat_to_unit([0.0, 5.0]),
        ];
        let (c, r) = enclosing_cap(pts.into_iter(), lonlat_to_unit([0.0, 4.0]));
        assert_abs_diff_eq!(r.to_degrees(), 10.0, epsilon = 0.05);
        assert_abs_diff_eq!(unit_to_lonlat(c)[1], 0.0, epsilon = 0.1);
    }

    #[test]
    fn contains_point_antipodal_to_reference_uses_two_legs() {
        // Big polygon (AsGiven complement): reference sits near the square,
        // so points on the far side are nearly antipodal to it.
        let comp = Polygon::from_lonlat_with(
            &[vec![[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0]]],
            Orientation::AsGiven,
        )
        .unwrap();
        assert!(comp.contains_lonlat([-179.5, -0.5]));
        assert!(comp.contains_lonlat([-179.5, -0.6]));
    }
}
