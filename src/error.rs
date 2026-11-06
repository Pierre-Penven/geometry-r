//! Crate-wide error type.

use thiserror::Error;

/// Errors returned by fallible operations of this crate.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum GeoError {
    /// A latitude outside `[-90, 90]` degrees, or a non-finite coordinate.
    #[error("invalid coordinate: {0}")]
    InvalidCoordinate(String),

    /// An argument outside its documented domain.
    #[error("invalid argument: {0}")]
    InvalidArgument(String),

    /// A polygon ring has fewer than three distinct vertices or is otherwise
    /// degenerate.
    #[error("degenerate polygon: {0}")]
    DegeneratePolygon(String),

    /// An iterative algorithm did not converge (e.g. Vincenty for nearly
    /// antipodal points).
    #[error("algorithm did not converge: {0}")]
    NoConvergence(&'static str),

    /// A ray (line of sight, beam edge, boresight) does not intersect the
    /// Earth.
    #[error("ray does not intersect the earth")]
    NoIntersection,

    /// The GeoJSON input is malformed or uses an unsupported geometry type.
    #[error("geojson: {0}")]
    GeoJson(String),
}

/// Convenience alias for `Result<T, GeoError>`.
pub type Result<T> = std::result::Result<T, GeoError>;
