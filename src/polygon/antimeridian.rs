//! Conversion of spherical polygons to RFC 7946 lon/lat polygons (split at
//! the antimeridian, closed through the poles) and antimeridian-aware
//! bounding boxes.

use super::planar::{orient_shape, push_window_offset, Contour};
use super::{lonlat_to_unit, ring_edges, unit_to_lonlat, LonLat, Polygon};
use crate::math::{wrap_lon_deg_closed, wrap_lon_delta_deg, Vec3};
use i_overlay::core::{fill_rule::FillRule, overlay_rule::OverlayRule};
use i_overlay::float::single::SingleFloatOverlay;

/// A lon/lat polygon: closed rings, exterior counter-clockwise first, holes
/// clockwise.
pub type LonLatPolygon = Vec<Vec<LonLat>>;

/// Sample point `[lon, lat]` used to fix the winding offset of the world
/// window (any point works; this one is unlikely to sit on an edge).
const WORLD_TEST_POINT: LonLat = [12.345_678_9, 7.123_456_7];

/// The world rectangle, clockwise in mirrored coordinates (filled).
const WORLD_FILLED: [[f64; 2]; 4] = [[-180.0, -90.0], [-180.0, 90.0], [180.0, 90.0], [180.0, -90.0]];

/// Output coordinates within this many degrees of ±180 / ±90 are snapped
/// to the exact value (planar snapping leaves ~1e-7° of noise).
const SNAP_DEG: f64 = 1e-7;

/// One ring after crossing insertion and longitude unwrapping.
struct Unwrapped {
    pts: Vec<LonLat>,
    /// Net number of turns around the polar axis (−1, 0 or +1).
    turns: i32,
    crosses: bool,
}

/// Inserts the exact great-circle crossing points with the antimeridian and
/// unwraps longitudes so that consecutive vertices differ by at most 180°.
/// Rings that turn around the pole are closed through it.
fn unwrap_ring(ring: &[Vec3]) -> Unwrapped {
    let mut pts: Vec<LonLat> = Vec::with_capacity(ring.len() + 4);
    let mut crosses = false;
    let push = |ll: LonLat, pts: &mut Vec<LonLat>| {
        let lon = match pts.last() {
            None => ll[0],
            Some(prev) => prev[0] + wrap_lon_delta_deg(ll[0] - prev[0]),
        };
        pts.push([lon, ll[1]]);
    };
    for (a, b) in ring_edges(ring) {
        push(unit_to_lonlat(a), &mut pts);
        if let Some(lat) = antimeridian_crossing_lat(a, b) {
            crosses = true;
            push([180.0, lat], &mut pts);
        }
    }
    let first = pts[0][0];
    let last = pts[pts.len() - 1][0];
    let end = last + wrap_lon_delta_deg(unit_to_lonlat(ring[0])[0] - last);
    let turns = ((end - first) / 360.0).round() as i32;
    if turns != 0 {
        let pole = if turns > 0 { 90.0 } else { -90.0 };
        let lat0 = pts[0][1];
        pts.extend([[end, lat0], [end, pole], [first, pole]]);
    }
    Unwrapped { pts, turns, crosses }
}

/// Latitude where the minor arc `a → b` crosses the antimeridian, if it does.
fn antimeridian_crossing_lat(a: Vec3, b: Vec3) -> Option<f64> {
    let straddles_y0 = (a.y > 0.0 && b.y < 0.0) || (a.y < 0.0 && b.y > 0.0);
    if !straddles_y0 {
        return None;
    }
    let p = (a * b.y.abs() + b * a.y.abs()).normalize();
    (p.x < 0.0).then(|| p.z.atan2(p.x.hypot(p.y)).to_degrees())
}

fn snap(v: f64, target: f64) -> f64 {
    if (v - target).abs() < SNAP_DEG {
        target
    } else {
        v
    }
}

/// Repeats the first vertex at the end (GeoJSON rings are closed).
fn close(mut r: Vec<LonLat>) -> Vec<LonLat> {
    r.push(r[0]);
    r
}

/// Converts a spherical polygon to one or more RFC 7946 lon/lat polygons.
///
/// * Polygons crossing the antimeridian are split into pieces that touch
///   ±180°.
/// * Rings that enclose a pole are closed along the pole latitude.
/// * Polygons that need no split are written exactly; the others go through
///   a planar clip whose snapping error is below 1e-6° (≈ 10 cm).
///
/// `max_edge_deg` densifies edges first so that great-circle edges render
/// faithfully in lon/lat viewers.
#[must_use]
pub fn to_rfc7946(poly: &Polygon, max_edge_deg: Option<f64>) -> Vec<LonLatPolygon> {
    let densified;
    let poly = match max_edge_deg.filter(|d| *d > 0.0) {
        Some(d) => {
            densified = poly
                .densify(d.to_radians())
                .expect("positive densification step is valid");
            &densified
        }
        None => poly,
    };
    let rings: Vec<Unwrapped> = poly.rings().iter().map(|r| unwrap_ring(r)).collect();
    let needs_split = !poly.is_small() || rings.iter().any(|r| r.crosses || r.turns != 0);
    if needs_split {
        split_with_overlay(poly, &rings)
    } else {
        vec![exact_rings(poly)]
    }
}

/// Rings written vertex for vertex (no split needed).
fn exact_rings(poly: &Polygon) -> LonLatPolygon {
    poly.rings()
        .iter()
        .map(|r| close(r.iter().map(|&v| unit_to_lonlat(v)).collect()))
        .collect()
}

/// Copies of every unwrapped ring shifted by multiples of 360° that overlap
/// the world window, in mirrored coordinates.
fn shifted_copies(rings: &[Unwrapped]) -> Vec<Contour> {
    let mut out = Vec::new();
    for r in rings {
        let (lo, hi) = r
            .pts
            .iter()
            .fold((f64::MAX, f64::MIN), |(lo, hi), p| (lo.min(p[0]), hi.max(p[0])));
        for k in -2..=2 {
            let shift = 360.0 * f64::from(k);
            if hi + shift > -180.0 && lo + shift < 180.0 {
                out.push(r.pts.iter().map(|p| [p[0] + shift, -p[1]]).collect());
            }
        }
    }
    out
}

fn split_with_overlay(poly: &Polygon, rings: &[Unwrapped]) -> Vec<LonLatPolygon> {
    let mut subject = shifted_copies(rings);
    let originals = subject.clone();
    let restorer = Restorer::new(&originals);
    let test_mirrored = [WORLD_TEST_POINT[0], -WORLD_TEST_POINT[1]];
    let want = poly.contains(lonlat_to_unit(WORLD_TEST_POINT));
    push_window_offset(&mut subject, &WORLD_FILLED, test_mirrored, want);
    let clip = vec![WORLD_FILLED.to_vec()];
    subject
        .overlay(&clip, OverlayRule::Intersect, FillRule::Positive)
        .into_iter()
        .map(|mut shape| {
            orient_shape(&mut shape);
            shape.iter().map(|c| close(unmirror(&restorer.restore(c)))).collect()
        })
        .collect()
}

/// Mirrored planar contour → lon/lat ring, snapping to ±180 / ±90.
fn unmirror(c: &[[f64; 2]]) -> Vec<LonLat> {
    c.iter()
        .map(|p| [snap(snap(p[0], 180.0), -180.0), snap(snap(-p[1], 90.0), -90.0)])
        .collect()
}

/// Re-inserts input vertices that the planar overlay dropped because they
/// were collinear in the lon/lat plane (e.g. points along one parallel).
/// They are not collinear on the sphere, so dropping them would change the
/// great-circle geometry.
struct Restorer<'a> {
    contours: &'a [Vec<[f64; 2]>],
    cells: std::collections::HashMap<(i64, i64), Vec<(usize, usize)>>,
}

impl<'a> Restorer<'a> {
    const CELL: f64 = 1e-5;

    fn key(p: [f64; 2]) -> (i64, i64) {
        ((p[0] / Self::CELL).round() as i64, (p[1] / Self::CELL).round() as i64)
    }

    fn new(contours: &'a [Vec<[f64; 2]>]) -> Self {
        let mut cells: std::collections::HashMap<(i64, i64), Vec<(usize, usize)>> = std::collections::HashMap::new();
        for (ci, c) in contours.iter().enumerate() {
            for (vi, &p) in c.iter().enumerate() {
                cells.entry(Self::key(p)).or_default().push((ci, vi));
            }
        }
        Self { contours, cells }
    }

    /// All input vertices within one cell of `p`.
    fn near(&self, p: [f64; 2]) -> impl Iterator<Item = (usize, usize)> + '_ {
        let (kx, ky) = Self::key(p);
        (-1..=1)
            .flat_map(move |dx| (-1..=1).map(move |dy| (kx + dx, ky + dy)))
            .filter_map(|k| self.cells.get(&k))
            .flatten()
            .copied()
            .filter(move |&(ci, vi)| {
                let q = self.contours[ci][vi];
                (q[0] - p[0]).hypot(q[1] - p[1]) < Self::CELL
            })
    }

    /// Parameter of `p` along `a → b`, if `p` lies strictly inside the segment.
    fn seg_param(a: [f64; 2], b: [f64; 2], p: [f64; 2]) -> Option<f64> {
        let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
        let l2 = dx * dx + dy * dy;
        if l2 == 0.0 {
            return None;
        }
        let t = ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / l2;
        let dist = ((p[0] - a[0]) * dy - (p[1] - a[1]) * dx).abs() / l2.sqrt();
        let end_gap = Self::CELL / l2.sqrt();
        (t > end_gap && t < 1.0 - end_gap && dist < Self::CELL).then_some(t)
    }

    /// Walks input contours from `a` towards `b`, collecting the vertices
    /// that lie on the segment (switching contour where copies touch).
    fn restore_edge(&self, a: [f64; 2], b: [f64; 2], out: &mut Vec<[f64; 2]>) {
        let mut cur = a;
        let mut t_cur = 0.0;
        while let Some((q, t)) = self.next_on_segment(a, b, cur, t_cur) {
            out.push(q);
            cur = q;
            t_cur = t;
        }
    }

    /// Closest input vertex beyond parameter `t_cur` on `a → b` among the
    /// contour neighbours of the vertices at `cur`.
    fn next_on_segment(&self, a: [f64; 2], b: [f64; 2], cur: [f64; 2], t_cur: f64) -> Option<([f64; 2], f64)> {
        self.near(cur)
            .flat_map(|(ci, vi)| {
                let c = &self.contours[ci];
                let len = c.len();
                [c[(vi + 1) % len], c[(vi + len - 1) % len]]
            })
            .filter_map(|q| Self::seg_param(a, b, q).map(|t| (q, t)))
            .filter(|&(_, t)| t > t_cur)
            .min_by(|x, y| x.1.total_cmp(&y.1))
    }

    fn restore(&self, contour: &[[f64; 2]]) -> Vec<[f64; 2]> {
        let n = contour.len();
        let mut out = Vec::with_capacity(n);
        for k in 0..n {
            let (a, b) = (contour[k], contour[(k + 1) % n]);
            out.push(a);
            self.restore_edge(a, b, &mut out);
        }
        out
    }
}

/// An RFC 7946 bounding box in degrees. When the box crosses the
/// antimeridian, `west > east`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LonLatBox {
    /// Western longitude.
    pub west: f64,
    /// Southern latitude.
    pub south: f64,
    /// Eastern longitude.
    pub east: f64,
    /// Northern latitude.
    pub north: f64,
}

impl LonLatBox {
    /// True when the box crosses the antimeridian.
    #[must_use]
    pub fn crosses_antimeridian(&self) -> bool {
        self.west > self.east
    }

    /// True when the `[lon, lat]` point is inside the box.
    #[must_use]
    pub fn contains(&self, ll: LonLat) -> bool {
        let lat_ok = ll[1] >= self.south && ll[1] <= self.north;
        let lon = wrap_lon_deg_closed(ll[0]);
        let lon_ok = if self.crosses_antimeridian() {
            lon >= self.west || lon <= self.east
        } else {
            lon >= self.west && lon <= self.east
        };
        lat_ok && lon_ok
    }
}

/// Latitude extremes reached along the minor arc `a → b`.
fn arc_lat_range(a: Vec3, b: Vec3) -> (f64, f64) {
    let la = unit_to_lonlat(a)[1];
    let lb = unit_to_lonlat(b)[1];
    let (mut lo, mut hi) = (la.min(lb), la.max(lb));
    if let Some(n) = a.cross(b).try_normalize() {
        for z in [Vec3::Z, -Vec3::Z] {
            if let Some(m) = (z - n * n.dot(z)).try_normalize() {
                // m lies on the arc when it is between a and b.
                if a.cross(m).dot(n) > 0.0 && m.cross(b).dot(n) > 0.0 {
                    let lat = unit_to_lonlat(m)[1];
                    lo = lo.min(lat);
                    hi = hi.max(lat);
                }
            }
        }
    }
    (lo, hi)
}

/// Antimeridian- and pole-aware bounding box of a polygon.
#[must_use]
pub fn bounding_box(poly: &Polygon) -> LonLatBox {
    let ext = poly.exterior();
    let (mut south, mut north) = (90.0f64, -90.0f64);
    for (a, b) in ring_edges(ext) {
        let (lo, hi) = arc_lat_range(a, b);
        south = south.min(lo);
        north = north.max(hi);
    }
    let has_n = poly.contains(Vec3::Z);
    let has_s = poly.contains(-Vec3::Z);
    if has_n {
        north = 90.0;
    }
    if has_s {
        south = -90.0;
    }
    if has_n || has_s {
        return LonLatBox {
            west: -180.0,
            south,
            east: 180.0,
            north,
        };
    }
    // No pole inside: the exterior does not turn around the pole, so its
    // unwrapped longitudes span less than 360°.
    let u = unwrap_ring(ext);
    let (lo, hi) = u
        .pts
        .iter()
        .fold((f64::MAX, f64::MIN), |(lo, hi), p| (lo.min(p[0]), hi.max(p[0])));
    let mut west = wrap_lon_deg_closed(lo);
    let east = wrap_lon_deg_closed(hi);
    if west == 180.0 {
        west = -180.0;
    }
    LonLatBox {
        west,
        south,
        east,
        north,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::polygon::planar::signed_area;
    use crate::polygon::Orientation;
    use approx::assert_abs_diff_eq;

    fn rect(lon0: f64, lat0: f64, lon1: f64, lat1: f64) -> Polygon {
        Polygon::from_lonlat(&[vec![[lon0, lat0], [lon1, lat0], [lon1, lat1], [lon0, lat1]]]).unwrap()
    }

    fn reparse(parts: &[LonLatPolygon]) -> crate::polygon::MultiPolygon {
        crate::polygon::MultiPolygon::new(parts.iter().map(|p| Polygon::from_lonlat(p).unwrap()).collect())
    }

    #[test]
    fn simple_polygon_is_exact() {
        let p = rect(10.0, 10.0, 20.0, 20.0);
        let out = to_rfc7946(&p, None);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0][0].len(), 5);
        assert_eq!(out[0][0][0], out[0][0][4]);
        assert!(signed_area(&out[0][0][..4]) > 0.0, "exterior must be CCW");
    }

    #[test]
    fn split_at_antimeridian() {
        let p = rect(170.0, -10.0, -170.0, 10.0);
        let out = to_rfc7946(&p, None);
        assert_eq!(out.len(), 2);
        for part in &out {
            let lons: Vec<f64> = part[0].iter().map(|q| q[0]).collect();
            assert!(lons.iter().all(|&l| l >= 170.0 - 1e-9 || l <= -170.0 + 1e-9));
            assert!(lons.iter().any(|&l| l.abs() == 180.0));
            assert!(signed_area(&part[0]) > 0.0);
        }
        let back = reparse(&out).dissolve();
        assert_abs_diff_eq!(back.area_sr(), p.area_sr(), epsilon = 1e-9);
    }

    #[test]
    fn polar_cap_is_closed_through_pole() {
        let p = Polygon::circle(Vec3::Z, 20f64.to_radians(), 72).unwrap();
        let out = to_rfc7946(&p, None);
        assert_eq!(out.len(), 1);
        let ext = &out[0][0];
        assert!(ext.iter().any(|q| q[1] == 90.0));
        assert!(ext.iter().any(|q| q[0] == 180.0) && ext.iter().any(|q| q[0] == -180.0));
        // Planar area ≈ 360 × 20 (minus the chord bulges).
        assert!(signed_area(ext) > 360.0 * 19.0);
        let back = reparse(&out);
        assert!(back.contains_lonlat([0.0, 89.0]));
        // South polar cap.
        let s = Polygon::circle(-Vec3::Z, 10f64.to_radians(), 72).unwrap();
        let out = to_rfc7946(&s, None);
        assert!(out[0][0].iter().any(|q| q[1] == -90.0));
    }

    #[test]
    fn larger_than_hemisphere_output() {
        let comp = Polygon::from_lonlat_with(
            &[vec![[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0]]],
            Orientation::AsGiven,
        )
        .unwrap();
        let out = to_rfc7946(&comp, None);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].len(), 2, "world rectangle with one hole");
    }

    #[test]
    fn densify_output() {
        let p = rect(0.0, 0.0, 10.0, 10.0);
        let out = to_rfc7946(&p, Some(1.0));
        assert!(out[0][0].len() > 40);
    }

    #[test]
    fn bounding_boxes() {
        let b = bounding_box(&rect(170.0, -10.0, -170.0, 10.0));
        assert!(b.crosses_antimeridian());
        assert_abs_diff_eq!(b.west, 170.0, epsilon = 1e-9);
        assert_abs_diff_eq!(b.east, -170.0, epsilon = 1e-9);
        // Great-circle edges bulge poleward.
        assert!(b.north > 10.0 && b.north < 10.5);
        assert!(b.contains([180.0, 0.0]));
        assert!(!b.contains([0.0, 0.0]));
        let b = bounding_box(&Polygon::circle(Vec3::Z, 0.2, 36).unwrap());
        assert_eq!((b.west, b.east, b.north), (-180.0, 180.0, 90.0));
        assert!(b.contains([33.0, 89.0]));
        let b = bounding_box(&rect(10.0, 10.0, 20.0, 20.0));
        assert!(!b.crosses_antimeridian());
        assert!(b.contains([15.0, 15.0]) && !b.contains([25.0, 15.0]));
    }

    #[test]
    fn bbox_of_polygon_starting_on_antimeridian() {
        let p = Polygon::from_lonlat(&[vec![[180.0, 0.0], [-170.0, 0.0], [-170.0, 5.0], [180.0, 5.0]]]).unwrap();
        let b = bounding_box(&p);
        assert_eq!(b.west, -180.0);
        assert_abs_diff_eq!(b.east, -170.0, epsilon = 1e-9);
        assert!(!b.crosses_antimeridian());
    }

    #[test]
    fn bbox_of_south_polar_cap_reaches_minus_90() {
        let b = bounding_box(&Polygon::circle(-Vec3::Z, 0.3, 36).unwrap());
        assert_eq!((b.west, b.east, b.south), (-180.0, 180.0, -90.0));
        assert!(b.north < -70.0);
    }

    #[test]
    fn bbox_with_equatorial_edge_has_no_bulge() {
        let b = bounding_box(&rect(0.0, 0.0, 10.0, 10.0));
        assert_eq!(b.south, 0.0);
        assert!(b.north > 10.0);
    }

    #[test]
    fn bbox_full_world_contains_everything() {
        let b = LonLatBox {
            west: -180.0,
            south: -90.0,
            east: 180.0,
            north: 90.0,
        };
        assert!(b.contains([179.9, 0.0]) && b.contains([-180.0, 0.0]));
        assert!(!b.crosses_antimeridian());
    }

    #[test]
    fn seg_param_rejects_degenerate_and_outside_points() {
        assert_eq!(Restorer::seg_param([1.0, 1.0], [1.0, 1.0], [1.0, 1.0]), None);
        assert_eq!(Restorer::seg_param([0.0, 0.0], [10.0, 0.0], [11.0, 0.0]), None);
        assert_eq!(Restorer::seg_param([0.0, 0.0], [10.0, 0.0], [5.0, 1.0]), None);
        assert_abs_diff_eq!(Restorer::seg_param([0.0, 0.0], [10.0, 0.0], [5.0, 0.0]).unwrap(), 0.5);
    }

    #[test]
    fn crossing_latitude_only_for_arcs_through_antimeridian() {
        let a = lonlat_to_unit([170.0, 10.0]);
        let b = lonlat_to_unit([-170.0, 10.0]);
        assert!(
            antimeridian_crossing_lat(a, b).unwrap() > 10.0,
            "great circle bulges poleward"
        );
        // Crossing the prime meridian is not an antimeridian crossing.
        assert!(antimeridian_crossing_lat(lonlat_to_unit([-10.0, 0.0]), lonlat_to_unit([10.0, 0.0])).is_none());
        // Same side of y = 0.
        assert!(antimeridian_crossing_lat(lonlat_to_unit([10.0, 0.0]), lonlat_to_unit([20.0, 0.0])).is_none());
    }

    #[test]
    fn hole_across_antimeridian_becomes_notch_in_both_pieces() {
        let p = Polygon::from_lonlat(&[
            vec![[160.0, -20.0], [-160.0, -20.0], [-160.0, 20.0], [160.0, 20.0]],
            vec![[175.0, -5.0], [-175.0, -5.0], [-175.0, 5.0], [175.0, 5.0]],
        ])
        .unwrap();
        let out = to_rfc7946(&p, None);
        assert_eq!(out.len(), 2);
        // Each half-hole touches the cut, so it becomes a notch in the
        // piece's exterior rather than a separate hole ring.
        assert!(out.iter().all(|part| part.len() == 1));
        let back = reparse(&out).dissolve();
        assert_abs_diff_eq!(back.area_sr(), p.area_sr(), epsilon = 1e-8);
        assert!(!back.contains_lonlat([180.0, 0.0]));
    }

    #[test]
    fn restorer_picks_nearest_candidate_when_contours_touch() {
        // Two contours share (0, 0); both offer a vertex on the segment.
        let contours = vec![
            vec![[0.0, 0.0], [1.0, 0.0], [0.0, 5.0]],
            vec![[0.0, 0.0], [2.0, 0.0], [0.0, -5.0]],
        ];
        let r = Restorer::new(&contours);
        let mut out = Vec::new();
        r.restore_edge([0.0, 0.0], [3.0, 0.0], &mut out);
        assert_eq!(out, vec![[1.0, 0.0]]);
    }
}
