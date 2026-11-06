//! # geometry_r
//!
//! High-performance geodesy and spherical geometry on a round Earth.
//!
//! * **Coordinates**: geodetic ↔ ECEF (closed form), local ENU frames.
//! * **Distances**: great-circle, Vincenty, bearings, cross-track.
//! * **Satellite geometry**: look angles, body frames, scan angles in several
//!   conventions, inverse scan to the ground, nadir/elevation relations,
//!   GSO-arc separation, polarization skew, line of sight.
//! * **Footprints**: coverage contours, circular/elliptical beam footprints,
//!   analytic point-in-beam tests.
//! * **Spherical polygons**: area, containment, intersection, union,
//!   difference, xor, coverage fraction — correct across the antimeridian and
//!   around the poles.
//! * **GeoJSON** (RFC 7946) input and output with antimeridian splitting.
//! * **Link geometry** and **Sun** helpers (path loss, Doppler, sun outage,
//!   eclipse).
//!
//! Units: radians and metres internally; degrees only at the boundaries
//! (constructors, GeoJSON). Field names carry their unit.
//!
//! ```
//! use geometry_r::prelude::*;
//!
//! let station = Geodetic::from_degrees(48.85, 2.35, 35.0).unwrap();
//! let sat = geo_satellite_position(13.0_f64.to_radians());
//! let look = look_angles(&station, sat);
//! assert!(look.elevation_rad > 0.0);
//! ```

pub mod batch;
pub mod coords;
pub mod ellipsoid;
pub mod error;
pub mod footprint;
pub mod geodesic;
pub mod geojson;
pub mod link;
pub mod look;
pub mod math;
pub mod polygon;
pub mod satellite;
pub mod sun;

// Compile and run the README examples as doctests.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct ReadmeDoctests;

/// Commonly used items, for `use geometry_r::prelude::*`.
pub mod prelude {
    pub use crate::coords::{Ecef, Geodetic, LocalFrame, LookAngles};
    pub use crate::ellipsoid::{Ellipsoid, WGS84};
    pub use crate::error::{GeoError, Result};
    pub use crate::footprint::{beam_footprint, coverage_contour, Beam, BeamShape};
    pub use crate::geodesic::{great_circle_distance, initial_bearing, vincenty_inverse};
    pub use crate::look::{is_visible, look_angles};
    pub use crate::math::{Mat3, Vec3};
    pub use crate::polygon::{Cap, LonLat, MultiPolygon, Orientation, Polygon};
    pub use crate::satellite::{geo_satellite_position, SatelliteFrame, ScanAngles};
}
