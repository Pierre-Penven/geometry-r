//! Boolean operations on spherical polygons.
//!
//! # Algorithm
//! 1. **Single gnomonic chart** when every operand fits in a cap smaller
//!    than a hemisphere and all of them fit in one cap of at most 85°.
//!    The gnomonic projection maps great circles to straight lines, so the
//!    planar overlay (i_overlay) is exact for great-circle edges.
//! 2. **Cube-face fallback** otherwise: each ring is clipped to each of the
//!    six cube faces on the sphere (Sutherland–Hodgman against four
//!    great-circle planes), projected onto the face, and the winding offset
//!    of the face is fixed with one point-in-polygon test, so faces fully
//!    covered by a large region are handled. Faces are overlaid separately.
//! 3. The face pieces are merged in one chart when the result fits in one;
//!    otherwise they are returned as separate polygons (correct area and
//!    membership, with seams along cube edges).

use super::planar::{Contour, orient_shape, push_window_offset};
use super::{Cap, MultiPolygon, Polygon, clean_ring, enclosing_cap, ring_signed_area};
use crate::math::Vec3;
use i_overlay::core::{fill_rule::FillRule, overlay_rule::OverlayRule};
use i_overlay::float::single::SingleFloatOverlay;

/// Largest cap radius handled by the single-chart path. Beyond ~85° the
/// gnomonic scale grows quickly and planar snapping loses precision.
const CHART_MAX_RAD: f64 = 85.0 * std::f64::consts::PI / 180.0;

/// Sample point (face coordinates) used to fix the winding offset
/// of a cube face. Any interior point works; an "irrational-looking" one is
/// unlikely to sit exactly on a polygon edge.
const FACE_TEST_POINT: [f64; 2] = [0.123_456_789, -0.071_234_567];

/// The face square `[-1, 1]²`, counter-clockwise (filled).
const FACE_SQUARE_FILLED: [[f64; 2]; 4] = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];

/// A boolean operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BooleanOp {
    /// Points in both operands.
    Intersection,
    /// Points in either operand.
    Union,
    /// Points in the first operand but not the second.
    Difference,
    /// Points in exactly one operand.
    Xor,
}

impl BooleanOp {
    fn rule(self) -> OverlayRule {
        match self {
            Self::Intersection => OverlayRule::Intersect,
            Self::Union => OverlayRule::Union,
            Self::Difference => OverlayRule::Difference,
            Self::Xor => OverlayRule::Xor,
        }
    }
}

/// A gnomonic chart centred on `c`. The basis `e1, e2 = c × e1` is
/// right-handed seen from outside, so orientation is preserved.
#[derive(Debug, Clone, Copy)]
struct Chart {
    c: Vec3,
    e1: Vec3,
    e2: Vec3,
}

impl Chart {
    fn new(c: Vec3) -> Self {
        let e1 = c.any_orthogonal();
        Self { c, e1, e2: c.cross(e1) }
    }

    #[inline]
    fn project(&self, p: Vec3) -> [f64; 2] {
        let d = p.dot(self.c);
        [p.dot(self.e1) / d, p.dot(self.e2) / d]
    }

    #[inline]
    fn unproject(&self, q: [f64; 2]) -> Vec3 {
        (self.c + self.e1 * q[0] + self.e2 * q[1]).normalize()
    }

    fn project_polygons(&self, polys: &[&Polygon]) -> Vec<Contour> {
        polys
            .iter()
            .flat_map(|p| p.rings().iter())
            .map(|r| r.iter().map(|&v| self.project(v)).collect())
            .collect()
    }
}

/// Converts one i_overlay output shape back to a spherical polygon.
/// Returns `None` when the exterior collapses after unprojection.
fn shape_to_polygon(chart: &Chart, mut shape: Vec<Contour>) -> Option<Polygon> {
    orient_shape(&mut shape);
    let mut rings = Vec::with_capacity(shape.len());
    for (k, contour) in shape.into_iter().enumerate() {
        let ring = clean_ring(contour.iter().map(|&q| chart.unproject(q)).collect());
        if ring.len() >= 3 {
            rings.push(ring);
        } else if k == 0 {
            return None;
        }
    }
    (ring_signed_area(&rings[0]).abs() > 0.0).then(|| Polygon::build(rings))
}

/// Runs i_overlay in one chart and converts the result to polygons.
fn overlay_in_chart(chart: &Chart, subject: Vec<Contour>, clip: Vec<Contour>, rule: OverlayRule) -> Vec<Polygon> {
    subject
        .overlay(&clip, rule, FillRule::Positive)
        .into_iter()
        .filter_map(|shape| shape_to_polygon(chart, shape))
        .collect()
}

/// The six cube faces, each a chart whose face is the square `[-1, 1]²`.
fn cube_faces() -> [Chart; 6] {
    [Vec3::X, -Vec3::X, Vec3::Y, -Vec3::Y, Vec3::Z, -Vec3::Z].map(Chart::new)
}

/// Angular radius of the cap circumscribing one cube face.
fn face_cap_radius() -> f64 {
    (1.0f64 / 3f64.sqrt()).acos() + 1e-9
}

/// Sutherland–Hodgman clip of a spherical ring against the half-space
/// `n·p ≥ 0`.
fn clip_ring_plane(ring: &[Vec3], n: Vec3) -> Vec<Vec3> {
    let len = ring.len();
    let mut out = Vec::with_capacity(len + 4);
    for i in 0..len {
        let cur = ring[i];
        let prev = ring[(i + len - 1) % len];
        let (dc, dp) = (n.dot(cur), n.dot(prev));
        let (cin, pin) = (dc >= 0.0, dp >= 0.0);
        if cin != pin {
            // Point of the minor arc prev→cur on the plane.
            out.push((prev * dc.abs() + cur * dp.abs()).normalize());
        }
        if cin {
            out.push(cur);
        }
    }
    out
}

/// Clips a ring to a cube face and projects it onto the face square.
fn clip_ring_to_face(ring: &[Vec3], face: &Chart) -> Option<Contour> {
    let planes = [face.c - face.e1, face.c + face.e1, face.c - face.e2, face.c + face.e2];
    let mut r = ring.to_vec();
    for &n in &planes {
        if r.is_empty() {
            return None;
        }
        r = clip_ring_plane(&r, n);
    }
    (r.len() >= 3).then(|| {
        r.iter()
            .map(|&v| {
                let q = face.project(v);
                [q[0].clamp(-1.0, 1.0), q[1].clamp(-1.0, 1.0)]
            })
            .collect()
    })
}

/// Projected contours of `poly` on `face`, with the winding offset fixed so
/// that faces fully inside (or inside a hole of) a large region are right.
fn face_contours(poly: &Polygon, face: &Chart) -> Vec<Contour> {
    let mut contours: Vec<Contour> = poly.rings().iter().filter_map(|r| clip_ring_to_face(r, face)).collect();
    let want = poly.contains(face.unproject(FACE_TEST_POINT));
    push_window_offset(&mut contours, &FACE_SQUARE_FILLED, FACE_TEST_POINT, want);
    contours
}

/// Cap holding every polygon when all are small and fit in one chart.
fn fits_one_chart(polys: &[&Polygon]) -> Option<Cap> {
    if polys.is_empty() || !polys.iter().all(|p| p.is_small()) {
        return None;
    }
    // Each small polygon lies inside any cap (< 90°) holding its exterior
    // ring, so a cap over all exterior vertices bounds every region.
    let start = polys
        .iter()
        .fold(Vec3::ZERO, |acc, p| acc + p.cap().center())
        .try_normalize()
        .unwrap_or_else(|| polys[0].cap().center());
    let pts = polys.iter().flat_map(|p| p.exterior().iter().copied());
    let (c, r) = enclosing_cap(pts, start);
    (r <= CHART_MAX_RAD).then(|| Cap::new_unchecked(c, r))
}

fn overlay_on_cube_faces(subject: &[&Polygon], clip: &[&Polygon], rule: OverlayRule) -> Vec<Polygon> {
    let mut out = Vec::new();
    for face in cube_faces() {
        let fcap = Cap::new_unchecked(face.c, face_cap_radius());
        let contours_of = |set: &[&Polygon]| -> Vec<Contour> {
            set.iter()
                .filter(|p| p.cap().intersects(&fcap))
                .flat_map(|p| face_contours(p, &face))
                .collect()
        };
        let s = contours_of(subject);
        let c = contours_of(clip);
        let nothing_to_do =
            (s.is_empty() && c.is_empty()) || (rule == OverlayRule::Intersect && (s.is_empty() || c.is_empty()));
        if !nothing_to_do {
            out.extend(overlay_in_chart(&face, s, c, rule));
        }
    }
    out
}

/// Merges cube-face pieces in one chart when they fit in one.
fn dissolve_seams(pieces: Vec<Polygon>) -> Vec<Polygon> {
    if pieces.len() < 2 {
        return pieces;
    }
    let refs: Vec<&Polygon> = pieces.iter().collect();
    match fits_one_chart(&refs) {
        Some(cap) => {
            let chart = Chart::new(cap.center());
            overlay_in_chart(&chart, chart.project_polygons(&refs), Vec::new(), OverlayRule::Subject)
        }
        None => pieces,
    }
}

/// Overlays `subject` and `clip` polygon sets with an i_overlay rule.
pub(crate) fn overlay_sets(subject: &[&Polygon], clip: &[&Polygon], rule: OverlayRule) -> MultiPolygon {
    let all: Vec<&Polygon> = subject.iter().chain(clip.iter()).copied().collect();
    let polygons = match fits_one_chart(&all) {
        Some(cap) => {
            let chart = Chart::new(cap.center());
            overlay_in_chart(
                &chart,
                chart.project_polygons(subject),
                chart.project_polygons(clip),
                rule,
            )
        }
        None => dissolve_seams(overlay_on_cube_faces(subject, clip, rule)),
    };
    MultiPolygon::new(polygons)
}

impl MultiPolygon {
    fn refs(&self) -> Vec<&Polygon> {
        self.polygons.iter().collect()
    }

    fn caps_disjoint(&self, other: &MultiPolygon) -> bool {
        match (self.cap(), other.cap()) {
            (Some(a), Some(b)) => !a.intersects(&b),
            _ => true,
        }
    }

    /// Applies a boolean operation.
    #[must_use]
    pub fn boolean(&self, other: &MultiPolygon, op: BooleanOp) -> MultiPolygon {
        if self.caps_disjoint(other) {
            // Fast paths: no overlay needed.
            return match op {
                BooleanOp::Intersection => MultiPolygon::default(),
                BooleanOp::Difference => self.clone(),
                BooleanOp::Union | BooleanOp::Xor => {
                    let mut v = self.polygons.clone();
                    v.extend(other.polygons.iter().cloned());
                    MultiPolygon::new(v)
                }
            };
        }
        overlay_sets(&self.refs(), &other.refs(), op.rule())
    }

    /// Intersection.
    #[must_use]
    pub fn intersection(&self, other: &MultiPolygon) -> MultiPolygon {
        self.boolean(other, BooleanOp::Intersection)
    }

    /// Union.
    #[must_use]
    pub fn union(&self, other: &MultiPolygon) -> MultiPolygon {
        self.boolean(other, BooleanOp::Union)
    }

    /// Difference `self − other`.
    #[must_use]
    pub fn difference(&self, other: &MultiPolygon) -> MultiPolygon {
        self.boolean(other, BooleanOp::Difference)
    }

    /// Symmetric difference.
    #[must_use]
    pub fn xor(&self, other: &MultiPolygon) -> MultiPolygon {
        self.boolean(other, BooleanOp::Xor)
    }

    /// Merges overlapping or touching members (e.g. pieces of a polygon that
    /// was split at the antimeridian).
    #[must_use]
    pub fn dissolve(&self) -> MultiPolygon {
        if self.polygons.len() < 2 {
            return self.clone();
        }
        overlay_sets(&self.refs(), &[], OverlayRule::Subject)
    }

    /// Union of many multipolygons in a single overlay.
    #[must_use]
    pub fn union_all(items: &[MultiPolygon]) -> MultiPolygon {
        let refs: Vec<&Polygon> = items.iter().flat_map(|m| m.polygons.iter()).collect();
        if refs.len() < 2 {
            return MultiPolygon::new(refs.into_iter().cloned().collect());
        }
        overlay_sets(&refs, &[], OverlayRule::Subject)
    }

    /// Area (steradians) of `self ∩ other`.
    #[must_use]
    pub fn overlap_area_sr(&self, other: &MultiPolygon) -> f64 {
        self.intersection(other).area_sr()
    }

    /// Fraction of `self`'s area covered by the union of `covers`, in
    /// `[0, 1]`. Returns 0 for an empty `self`.
    #[must_use]
    pub fn coverage_fraction(&self, covers: &[MultiPolygon]) -> f64 {
        let total = self.area_sr();
        if total <= 0.0 {
            return 0.0;
        }
        let cover = MultiPolygon::union_all(covers);
        (self.intersection(&cover).area_sr() / total).clamp(0.0, 1.0)
    }
}

impl Polygon {
    /// Intersection with another polygon.
    #[must_use]
    pub fn intersection(&self, other: &Polygon) -> MultiPolygon {
        MultiPolygon::from(self.clone()).intersection(&other.clone().into())
    }

    /// Union with another polygon.
    #[must_use]
    pub fn union(&self, other: &Polygon) -> MultiPolygon {
        MultiPolygon::from(self.clone()).union(&other.clone().into())
    }

    /// Difference `self − other`.
    #[must_use]
    pub fn difference(&self, other: &Polygon) -> MultiPolygon {
        MultiPolygon::from(self.clone()).difference(&other.clone().into())
    }

    /// Symmetric difference.
    #[must_use]
    pub fn xor(&self, other: &Polygon) -> MultiPolygon {
        MultiPolygon::from(self.clone()).xor(&other.clone().into())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::polygon::{Orientation, lonlat_to_unit};
    use approx::assert_abs_diff_eq;
    use std::f64::consts::{FRAC_PI_2, PI};

    fn rect(lon0: f64, lat0: f64, lon1: f64, lat1: f64) -> Polygon {
        Polygon::from_lonlat(&[vec![[lon0, lat0], [lon1, lat0], [lon1, lat1], [lon0, lat1]]]).unwrap()
    }

    fn mp(p: Polygon) -> MultiPolygon {
        p.into()
    }

    fn check_identities(a: &MultiPolygon, b: &MultiPolygon, tol: f64) {
        let i = a.intersection(b).area_sr();
        let u = a.union(b).area_sr();
        let d = a.difference(b).area_sr();
        let x = a.xor(b).area_sr();
        let (aa, bb) = (a.area_sr(), b.area_sr());
        assert_abs_diff_eq!(u, aa + bb - i, epsilon = tol);
        assert_abs_diff_eq!(d, aa - i, epsilon = tol);
        assert_abs_diff_eq!(x, u - i, epsilon = tol);
    }

    #[test]
    fn simple_overlap() {
        let a = mp(rect(0.0, 0.0, 2.0, 2.0));
        let b = mp(rect(1.0, 1.0, 3.0, 3.0));
        let i = a.intersection(&b);
        assert_eq!(i.polygons.len(), 1);
        assert!(i.contains_lonlat([1.5, 1.5]));
        assert!(!i.contains_lonlat([0.5, 0.5]));
        check_identities(&a, &b, 1e-10);
        // Intersection of a great-circle square with itself is itself.
        assert_abs_diff_eq!(a.intersection(&a).area_sr(), a.area_sr(), epsilon = 1e-12);
    }

    #[test]
    fn across_antimeridian() {
        let a = mp(rect(170.0, -5.0, -170.0, 5.0));
        let b = mp(rect(175.0, -10.0, -175.0, 0.0));
        let i = a.intersection(&b);
        assert_eq!(i.polygons.len(), 1);
        assert!(i.contains_lonlat([180.0, -2.0]));
        check_identities(&a, &b, 1e-10);
        // Pieces split at ±180 dissolve into one polygon.
        let west = rect(170.0, -5.0, 180.0, 5.0);
        let east = rect(-180.0, -5.0, -170.0, 5.0);
        let parts = west.area_sr() + east.area_sr();
        let d = MultiPolygon::new(vec![west, east]).dissolve();
        assert_eq!(d.polygons.len(), 1);
        assert_abs_diff_eq!(d.area_sr(), parts, epsilon = 1e-9);
        // The shared meridian edge disappears; the 6 corners remain.
        assert_eq!(d.polygons[0].exterior().len(), 6);
    }

    #[test]
    fn around_the_pole() {
        let cap = Polygon::circle(Vec3::Z, 15f64.to_radians(), 360).unwrap();
        let wedge = Polygon::from_lonlat(&[vec![[0.0, 70.0], [90.0, 70.0], [0.0, 89.0]]]).unwrap();
        let a = mp(cap);
        let b = mp(wedge);
        check_identities(&a, &b, 1e-10);
        let d = a.difference(&b);
        assert!(d.contains_lonlat([180.0, 85.0]));
        assert!(!d.contains_lonlat([30.0, 80.0]));
    }

    #[test]
    fn holes_are_created_and_kept() {
        let a = mp(rect(0.0, 0.0, 10.0, 10.0));
        let b = mp(rect(3.0, 3.0, 6.0, 6.0));
        let d = a.difference(&b);
        assert_eq!(d.polygons.len(), 1);
        assert_eq!(d.polygons[0].holes().len(), 1);
        assert!(!d.contains_lonlat([4.0, 4.0]));
        assert!(d.contains_lonlat([1.0, 1.0]));
        let u = d.union(&b);
        assert_abs_diff_eq!(u.area_sr(), a.area_sr(), epsilon = 1e-10);
    }

    #[test]
    fn disjoint_fast_paths() {
        let a = mp(rect(0.0, 0.0, 1.0, 1.0));
        let b = mp(rect(50.0, 0.0, 51.0, 1.0));
        assert!(a.intersection(&b).is_empty());
        assert_eq!(a.union(&b).polygons.len(), 2);
        assert_eq!(a.difference(&b), a);
        assert_eq!(a.xor(&b).polygons.len(), 2);
        let e = MultiPolygon::default();
        assert!(e.intersection(&a).is_empty());
        assert_eq!(a.union(&e), a);
    }

    #[test]
    fn large_regions_use_cube_faces() {
        // Two GEO-like coverage caps (radius 81°) 60° apart.
        let r = 81f64.to_radians();
        let a = mp(Polygon::circle(lonlat_to_unit([0.0, 0.0]), r, 720).unwrap());
        let b = mp(Polygon::circle(lonlat_to_unit([60.0, 0.0]), r, 720).unwrap());
        check_identities(&a, &b, 1e-6);
        let i = a.intersection(&b);
        assert!(i.contains_lonlat([30.0, 0.0]));
        assert!(!i.contains_lonlat([-70.0, 0.0]));
        let u = a.union(&b);
        assert!(u.contains_lonlat([-70.0, 0.0]));
        assert!(u.contains_lonlat([130.0, 0.0]));
        assert!(!u.contains_lonlat([-150.0, 0.0]));
        // Lens area by numeric integration of the two caps is hard; use symmetry:
        // intersection is symmetric about lon 30.
        let ii = b.intersection(&a);
        assert_abs_diff_eq!(i.area_sr(), ii.area_sr(), epsilon = 1e-6);
    }

    #[test]
    fn larger_than_hemisphere_operands() {
        // Complement of a small square (AsGiven) ∩ a big cap.
        let comp = Polygon::from_lonlat_with(
            &[vec![[0.0, 0.0], [0.0, 2.0], [2.0, 2.0], [2.0, 0.0]]],
            Orientation::AsGiven,
        )
        .unwrap();
        let sq = rect(0.0, 0.0, 2.0, 2.0);
        assert_abs_diff_eq!(comp.area_sr() + sq.area_sr(), 4.0 * PI, epsilon = 1e-12);
        let disk = mp(Polygon::circle(lonlat_to_unit([1.0, 1.0]), 0.5, 360).unwrap());
        let r = disk.intersection(&mp(comp.clone()));
        assert_abs_diff_eq!(r.area_sr(), disk.area_sr() - sq.area_sr(), epsilon = 1e-9);
        assert!(!r.contains_lonlat([1.0, 1.0]));
        assert!(r.contains_lonlat([1.0, 10.0]));
        // Complement ∪ square = whole sphere area.
        let u = mp(comp.clone()).union(&mp(sq.clone()));
        assert_abs_diff_eq!(u.area_sr(), 4.0 * PI, epsilon = 1e-6);
        // Complement − disk.
        let d = mp(comp).difference(&disk);
        assert_abs_diff_eq!(d.area_sr(), 4.0 * PI - disk.area_sr(), epsilon = 1e-6);
    }

    #[test]
    fn hemisphere_scale_intersection_matches_formula() {
        // Two hemispheres-ish caps of 89° centred 90° apart: lune-like region.
        let a = mp(Polygon::circle(Vec3::Z, FRAC_PI_2 - 0.01, 2000).unwrap());
        let b = mp(Polygon::circle(Vec3::X, FRAC_PI_2 - 0.01, 2000).unwrap());
        let i = a.intersection(&b).area_sr();
        // Two hemispheres at 90° meet in a lune of area π; shrinking each cap
        // by 0.01 rad removes a strip of width 0.01 along both half circles.
        assert_abs_diff_eq!(i, PI - 2.0 * PI * 0.01, epsilon = 1e-3);
        check_identities(&a, &b, 1e-6);
    }

    #[test]
    fn union_all_and_coverage() {
        let target = mp(rect(0.0, 0.0, 10.0, 10.0));
        let beams: Vec<MultiPolygon> = vec![mp(rect(0.0, 0.0, 5.0, 10.0)), mp(rect(4.0, 0.0, 10.0, 5.0))];
        let f = target.coverage_fraction(&beams);
        // Covered: left half + bottom-right quarter ≈ 75 %.
        assert_abs_diff_eq!(f, 0.75, epsilon = 0.01);
        assert_eq!(MultiPolygon::default().coverage_fraction(&beams), 0.0);
        assert_eq!(target.coverage_fraction(&[]), 0.0);
        let u = MultiPolygon::union_all(&beams);
        assert_eq!(u.polygons.len(), 1);
        assert_abs_diff_eq!(target.overlap_area_sr(&u), u.area_sr(), epsilon = 1e-9);
        assert_eq!(MultiPolygon::union_all(&beams[..1]), beams[0]);
    }

    #[test]
    fn polygon_convenience_methods() {
        let a = rect(0.0, 0.0, 2.0, 2.0);
        let b = rect(1.0, 0.0, 3.0, 2.0);
        assert_abs_diff_eq!(a.intersection(&b).area_sr() * 2.0, a.area_sr(), epsilon = 1e-6);
        assert!(a.union(&b).area_sr() > a.area_sr());
        assert!(a.difference(&b).area_sr() < a.area_sr());
        assert!(a.xor(&b).area_sr() > 0.0);
        let one = MultiPolygon::from(a.clone());
        assert_eq!(one.dissolve(), one);
    }

    #[test]
    fn cube_path_skips_faces_with_nothing_on_them() {
        // A 100° cap forces the cube path; the face opposite to it is empty.
        let big = mp(Polygon::circle(lonlat_to_unit([0.0, 0.0]), 100f64.to_radians(), 720).unwrap());
        let small = mp(rect(-1.0, -1.0, 1.0, 1.0));
        let u = big.union(&small);
        assert_abs_diff_eq!(u.area_sr(), big.area_sr(), epsilon = 1e-6);
        assert!(!u.contains_lonlat([180.0, 0.0]));
        let d = big.difference(&small);
        assert_abs_diff_eq!(d.area_sr(), big.area_sr() - small.area_sr(), epsilon = 1e-6);
    }

    #[test]
    fn shape_to_polygon_drops_collapsed_exterior() {
        let chart = Chart::new(Vec3::X);
        let tiny = vec![vec![[0.0, 0.0], [0.0, 1e-14], [1e-14, 1e-14], [1e-14, 0.0]]];
        assert!(shape_to_polygon(&chart, tiny).is_none());
    }

    #[test]
    fn shape_to_polygon_drops_collapsed_hole_but_keeps_exterior() {
        let chart = Chart::new(Vec3::X);
        let shape = vec![
            vec![[-0.1, -0.1], [-0.1, 0.1], [0.1, 0.1], [0.1, -0.1]],
            vec![[0.0, 0.0], [1e-14, 0.0], [1e-14, 1e-14], [0.0, 1e-14]],
        ];
        let p = shape_to_polygon(&chart, shape).unwrap();
        assert!(p.holes().is_empty());
        assert!(p.contains(Vec3::X));
    }

    #[test]
    fn chart_round_trip_and_orientation() {
        let chart = Chart::new(lonlat_to_unit([30.0, 40.0]));
        let v = lonlat_to_unit([31.0, 41.0]);
        let q = chart.project(v);
        assert_abs_diff_eq!(chart.unproject(q).dot(v), 1.0, epsilon = 1e-14);
        // Counter-clockwise on the sphere stays counter-clockwise (filled).
        let ring: Vec<[f64; 2]> = rect(30.0, 40.0, 31.0, 41.0)
            .exterior()
            .iter()
            .map(|&p| chart.project(p))
            .collect();
        assert!(crate::polygon::planar::signed_area(&ring) > 0.0);
    }

    #[test]
    fn clip_ring_to_face_misses_far_ring() {
        let face = Chart::new(Vec3::X);
        let far = rect(170.0, -1.0, 172.0, 1.0);
        assert!(clip_ring_to_face(far.exterior(), &face).is_none());
        let near = rect(-1.0, -1.0, 1.0, 1.0);
        let c = clip_ring_to_face(near.exterior(), &face).unwrap();
        assert!(c.iter().all(|q| q[0].abs() <= 1.0 && q[1].abs() <= 1.0));
    }
}
