//! Radio-link geometry: slant range, delay, free-space path loss, Doppler and
//! power flux density.

use crate::coords::Ecef;
use crate::ellipsoid::SPEED_OF_LIGHT_M_S;
use crate::math::Vec3;
use std::f64::consts::PI;

/// Straight-line distance between two ECEF points, in metres.
#[inline]
#[must_use]
pub fn slant_range(a: Ecef, b: Ecef) -> f64 {
    a.distance(b)
}

/// One-way propagation delay in seconds over `range_m` metres (vacuum).
#[inline]
#[must_use]
pub fn propagation_delay_s(range_m: f64) -> f64 {
    range_m / SPEED_OF_LIGHT_M_S
}

/// Free-space path loss in dB: `20·log10(4π d f / c)`.
#[inline]
#[must_use]
pub fn free_space_path_loss_db(range_m: f64, frequency_hz: f64) -> f64 {
    20.0 * (4.0 * PI * range_m * frequency_hz / SPEED_OF_LIGHT_M_S).log10()
}

/// Power flux density at `range_m` in dBW/m² for a given EIRP (dBW):
/// `EIRP − 10·log10(4π d²)`.
#[inline]
#[must_use]
pub fn power_flux_density_dbw_m2(eirp_dbw: f64, range_m: f64) -> f64 {
    eirp_dbw - 10.0 * (4.0 * PI * range_m * range_m).log10()
}

/// Range rate (m/s, positive when the distance increases) between an
/// observer and a target, both with ECEF positions and velocities expressed
/// in the same frame.
#[must_use]
pub fn range_rate(obs_pos: Ecef, obs_vel: Vec3, target_pos: Ecef, target_vel: Vec3) -> f64 {
    let d = target_pos - obs_pos;
    let n = d.norm();
    if n == 0.0 {
        return 0.0;
    }
    d.dot(target_vel - obs_vel) / n
}

/// Doppler shift in Hz of a carrier at `frequency_hz` for a given range rate
/// (non-relativistic): `−f · ṙ / c`. Approaching targets give positive shifts.
#[inline]
#[must_use]
pub fn doppler_shift_hz(frequency_hz: f64, range_rate_m_s: f64) -> f64 {
    -frequency_hz * range_rate_m_s / SPEED_OF_LIGHT_M_S
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn path_loss_and_delay() {
        // GEO at 12 GHz: ~205.1 dB.
        assert_abs_diff_eq!(free_space_path_loss_db(35_786e3, 12e9), 205.1, epsilon = 0.1);
        assert_abs_diff_eq!(propagation_delay_s(SPEED_OF_LIGHT_M_S), 1.0);
        assert_abs_diff_eq!(slant_range(Vec3::ZERO, Vec3::new(3.0, 4.0, 0.0)), 5.0);
    }

    #[test]
    fn pfd() {
        // 50 dBW EIRP at GEO range: ~-112 dBW/m².
        assert_abs_diff_eq!(power_flux_density_dbw_m2(50.0, 35_786e3), -112.07, epsilon = 0.05);
    }

    #[test]
    fn doppler() {
        let rr = range_rate(
            Vec3::ZERO,
            Vec3::ZERO,
            Vec3::new(1000.0, 0.0, 0.0),
            Vec3::new(-7000.0, 0.0, 0.0),
        );
        assert_abs_diff_eq!(rr, -7000.0);
        assert!(doppler_shift_hz(2e9, rr) > 0.0);
        assert_abs_diff_eq!(doppler_shift_hz(2e9, rr), 46_699.0, epsilon = 0.1);
        assert_eq!(range_rate(Vec3::X, Vec3::ZERO, Vec3::X, Vec3::Y), 0.0);
        // Purely tangential motion gives no range rate.
        assert_abs_diff_eq!(range_rate(Vec3::ZERO, Vec3::ZERO, Vec3::X, Vec3::Y), 0.0);
    }
}
