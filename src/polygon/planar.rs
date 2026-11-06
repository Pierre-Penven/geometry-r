//! Planar helpers shared by the chart overlay (`boolean`) and the RFC 7946
//! splitter (`antimeridian`): signed area, winding number, the winding-offset
//! correction, and the adapter that hides i_overlay's orientation convention.
//!
//! Convention: a contour that encloses "inside" is counter-clockwise, which
//! is what i_overlay ≥ 4 fills under `FillRule::Positive`. Gnomonic charts
//! preserve orientation, so a spherical exterior ring (counter-clockwise seen
//! from outside) projects to a filled contour without any flipping.

/// An open planar contour (first point not repeated).
pub(crate) type Contour = Vec<[f64; 2]>;

/// Shoelace signed area, positive for counter-clockwise contours.
pub(crate) fn signed_area(c: &[[f64; 2]]) -> f64 {
    let n = c.len();
    (0..n)
        .map(|i| {
            let (a, b) = (c[i], c[(i + 1) % n]);
            a[0] * b[1] - b[0] * a[1]
        })
        .sum::<f64>()
        * 0.5
}

/// Winding number of a closed contour around `p` (counter-clockwise
/// positive), Sunday's crossing algorithm.
pub(crate) fn winding_number(c: &[[f64; 2]], p: [f64; 2]) -> i32 {
    let n = c.len();
    let mut w = 0;
    for i in 0..n {
        let (a, b) = (c[i], c[(i + 1) % n]);
        let cross = (b[0] - a[0]) * (p[1] - a[1]) - (p[0] - a[0]) * (b[1] - a[1]);
        if a[1] <= p[1] {
            if b[1] > p[1] && cross > 0.0 {
                w += 1;
            }
        } else if b[1] <= p[1] && cross < 0.0 {
            w -= 1;
        }
    }
    w
}

/// Fill count at `p`: each counter-clockwise (filled) contour around `p`
/// counts +1, each clockwise one −1.
pub(crate) fn fill_count(contours: &[Contour], p: [f64; 2]) -> i32 {
    contours.iter().map(|c| winding_number(c, p)).sum::<i32>()
}

/// Makes the fill count right everywhere in a convex window.
///
/// Clipping a ring to a window keeps the boundary pieces inside the window,
/// but loses the information of whether the window as a whole lies inside a
/// region whose boundary never enters it. The fill count is therefore off by
/// a constant `k` over the window. Comparing with one exact membership test
/// (`want`, 0 or 1) at a sample point gives `k`, and `|k|` copies of the window
/// contour (filled when `k > 0`, reversed when `k < 0`) fix it.
pub(crate) fn push_window_offset(contours: &mut Vec<Contour>, window_filled: &[[f64; 2]], p: [f64; 2], want: bool) {
    let k = i32::from(want) - fill_count(contours, p);
    for _ in 0..k.unsigned_abs() {
        if k > 0 {
            contours.push(window_filled.to_vec());
        } else {
            contours.push(window_filled.iter().rev().copied().collect());
        }
    }
}

/// Orients one i_overlay output shape: outer contour counter-clockwise,
/// holes clockwise. i_overlay 4 already does this (1.x did the opposite);
/// the check keeps us independent of that convention.
pub(crate) fn orient_shape(shape: &mut [Contour]) {
    for (k, contour) in shape.iter_mut().enumerate() {
        let a = signed_area(contour);
        if (k == 0 && a < 0.0) || (k > 0 && a > 0.0) {
            contour.reverse();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SQUARE_CCW: [[f64; 2]; 4] = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];

    fn cw() -> Contour {
        SQUARE_CCW.iter().rev().copied().collect()
    }

    fn ccw() -> Contour {
        SQUARE_CCW.to_vec()
    }

    #[test]
    fn signed_area_is_positive_for_ccw_and_negative_for_cw() {
        assert_eq!(signed_area(&SQUARE_CCW), 1.0);
        assert_eq!(signed_area(&cw()), -1.0);
    }

    #[test]
    fn winding_number_counts_orientation_inside_and_zero_outside() {
        assert_eq!(winding_number(&SQUARE_CCW, [0.5, 0.5]), 1);
        assert_eq!(winding_number(&cw(), [0.5, 0.5]), -1);
        assert_eq!(winding_number(&SQUARE_CCW, [1.5, 0.5]), 0);
        assert_eq!(winding_number(&SQUARE_CCW, [0.5, -0.5]), 0);
    }

    #[test]
    fn fill_count_treats_counter_clockwise_contours_as_filled() {
        assert_eq!(fill_count(&[ccw()], [0.5, 0.5]), 1);
        assert_eq!(fill_count(&[ccw(), ccw()], [0.5, 0.5]), 2);
        assert_eq!(fill_count(&[cw()], [0.5, 0.5]), -1);
        assert_eq!(fill_count(&[], [0.5, 0.5]), 0);
    }

    #[test]
    fn window_offset_adds_filled_window_when_point_should_be_inside() {
        let mut c: Vec<Contour> = Vec::new();
        push_window_offset(&mut c, &ccw(), [0.5, 0.5], true);
        assert_eq!(c.len(), 1);
        assert_eq!(fill_count(&c, [0.5, 0.5]), 1);
    }

    #[test]
    fn window_offset_adds_reversed_window_when_count_is_too_high() {
        let mut c = vec![ccw(), ccw()];
        push_window_offset(&mut c, &ccw(), [0.5, 0.5], true);
        assert_eq!(c.len(), 3);
        assert_eq!(fill_count(&c, [0.5, 0.5]), 1);
    }

    #[test]
    fn window_offset_adds_nothing_when_count_is_already_right() {
        let mut c = vec![ccw()];
        push_window_offset(&mut c, &ccw(), [0.5, 0.5], true);
        assert_eq!(c.len(), 1);
        let mut empty: Vec<Contour> = Vec::new();
        push_window_offset(&mut empty, &ccw(), [0.5, 0.5], false);
        assert!(empty.is_empty());
    }

    #[test]
    fn orient_shape_fixes_outer_and_hole_orientation() {
        let hole_small: Contour = vec![[0.2, 0.2], [0.4, 0.2], [0.4, 0.4], [0.2, 0.4]]; // counter-clockwise
        let mut shape = vec![cw(), hole_small];
        orient_shape(&mut shape);
        assert!(signed_area(&shape[0]) > 0.0, "outer must be counter-clockwise");
        assert!(signed_area(&shape[1]) < 0.0, "hole must be clockwise");
        // Already correct input is left unchanged.
        let before = shape.clone();
        orient_shape(&mut shape);
        assert_eq!(shape, before);
    }
}
