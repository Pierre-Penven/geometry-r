//! GeoJSON (RFC 7946) input and output for polygonal geometries.
//!
//! * Input: `Polygon`, `MultiPolygon`, `Feature`, `FeatureCollection` and
//!   `GeometryCollection` with polygonal members.
//! * Output: `Polygon` or `MultiPolygon`, split at the antimeridian and
//!   closed through the poles as RFC 7946 requires.

use crate::error::{GeoError, Result};
use crate::math::wrap_lon_delta_deg;
use crate::polygon::antimeridian::{to_rfc7946, LonLatPolygon};
use crate::polygon::{LonLat, MultiPolygon, Polygon};
use serde_json::{json, Map, Value};

/// How the straight edges of input GeoJSON rings are interpreted.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EdgeInterpretation {
    /// Edges are great-circle arcs (common for satellite-derived data).
    GreatCircle,
    /// Edges are straight lines in lon/lat (RFC 7946). They are densified so
    /// that no step exceeds `max_step_deg` before conversion to great-circle
    /// edges.
    LonLatLinear {
        /// Maximum densification step in degrees.
        max_step_deg: f64,
    },
}

impl Default for EdgeInterpretation {
    fn default() -> Self {
        Self::LonLatLinear { max_step_deg: 1.0 }
    }
}

/// Options for reading GeoJSON.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct ReadOptions {
    /// Edge interpretation.
    pub edges: EdgeInterpretation,
    /// Merge overlapping or adjacent members (e.g. antimeridian-split
    /// pieces) after reading.
    pub dissolve: bool,
}

/// Options for writing GeoJSON.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct WriteOptions {
    /// Densify great-circle edges so that no edge exceeds this many degrees,
    /// for faithful rendering in lon/lat viewers.
    pub densify_deg: Option<f64>,
}

fn err(msg: impl Into<String>) -> GeoError {
    GeoError::GeoJson(msg.into())
}

/// Parses a GeoJSON string into a [`MultiPolygon`].
///
/// # Errors
/// [`GeoError::GeoJson`] for malformed JSON, unsupported or non-polygonal
/// geometries; polygon construction errors otherwise.
pub fn parse_str(s: &str, opts: &ReadOptions) -> Result<MultiPolygon> {
    let v: Value = serde_json::from_str(s).map_err(|e| err(e.to_string()))?;
    parse_value(&v, opts)
}

/// Parses a GeoJSON value into a [`MultiPolygon`].
///
/// # Errors
/// See [`parse_str`].
pub fn parse_value(v: &Value, opts: &ReadOptions) -> Result<MultiPolygon> {
    let mut polys = Vec::new();
    collect(v, opts, &mut polys)?;
    let mp = MultiPolygon::new(polys);
    Ok(if opts.dissolve { mp.dissolve() } else { mp })
}

fn collect(v: &Value, opts: &ReadOptions, out: &mut Vec<Polygon>) -> Result<()> {
    let ty = v
        .get("type")
        .and_then(Value::as_str)
        .ok_or_else(|| err("missing \"type\""))?;
    match ty {
        "Polygon" => {
            out.push(parse_polygon(coords(v)?, opts)?);
        }
        "MultiPolygon" => {
            let arr = coords(v)?
                .as_array()
                .ok_or_else(|| err("MultiPolygon coordinates must be an array"))?;
            for p in arr {
                out.push(parse_polygon(p, opts)?);
            }
        }
        "Feature" => match v.get("geometry") {
            None | Some(Value::Null) => {}
            Some(g) => collect(g, opts, out)?,
        },
        "FeatureCollection" => {
            let feats = v
                .get("features")
                .and_then(Value::as_array)
                .ok_or_else(|| err("FeatureCollection without \"features\" array"))?;
            for f in feats {
                collect(f, opts, out)?;
            }
        }
        "GeometryCollection" => {
            let geoms = v
                .get("geometries")
                .and_then(Value::as_array)
                .ok_or_else(|| err("GeometryCollection without \"geometries\" array"))?;
            for g in geoms {
                collect(g, opts, out)?;
            }
        }
        other => {
            return Err(err(format!(
                "unsupported geometry type {other:?} (only polygonal types)"
            )))
        }
    }
    Ok(())
}

fn coords(v: &Value) -> Result<&Value> {
    v.get("coordinates").ok_or_else(|| err("missing \"coordinates\""))
}

fn parse_position(p: &Value) -> Result<LonLat> {
    let a = p.as_array().ok_or_else(|| err("position must be an array"))?;
    if a.len() < 2 {
        return Err(err("position needs at least 2 numbers"));
    }
    let lon = a[0].as_f64().ok_or_else(|| err("longitude must be a number"))?;
    let lat = a[1].as_f64().ok_or_else(|| err("latitude must be a number"))?;
    Ok([lon, lat])
}

fn parse_polygon(v: &Value, opts: &ReadOptions) -> Result<Polygon> {
    let rings = v
        .as_array()
        .ok_or_else(|| err("polygon coordinates must be an array of rings"))?;
    if rings.is_empty() {
        return Err(err("polygon without rings"));
    }
    let mut out = Vec::with_capacity(rings.len());
    for r in rings {
        let pts = r.as_array().ok_or_else(|| err("ring must be an array of positions"))?;
        let ring: Vec<LonLat> = pts.iter().map(parse_position).collect::<Result<_>>()?;
        out.push(match opts.edges {
            EdgeInterpretation::GreatCircle => ring,
            EdgeInterpretation::LonLatLinear { max_step_deg } => densify_lonlat(&ring, max_step_deg)?,
        });
    }
    Polygon::from_lonlat(&out)
}

/// Densifies a lon/lat ring with straight lon/lat steps.
///
/// Longitude jumps above 180° are read as crossing the antimeridian; an edge
/// from exactly −180° to +180° (or back) at constant latitude is a full
/// parallel, as used by RFC 7946 polar caps.
fn densify_lonlat(ring: &[LonLat], step: f64) -> Result<Vec<LonLat>> {
    if step.is_nan() || step <= 0.0 {
        return Err(GeoError::InvalidArgument("densification step must be > 0".into()));
    }
    let n = ring.len();
    let mut out = Vec::with_capacity(n * 2);
    for i in 0..n {
        let a = ring[i];
        let b = ring[(i + 1) % n];
        let raw = b[0] - a[0];
        let full_parallel = raw.abs() == 360.0 && a[1] == b[1];
        let dlon = if full_parallel { raw } else { wrap_lon_delta_deg(raw) };
        let dlat = b[1] - a[1];
        let k = (dlon.abs().max(dlat.abs()) / step).ceil().max(1.0) as usize;
        for j in 0..k {
            let t = j as f64 / k as f64;
            out.push([a[0] + dlon * t, a[1] + dlat * t]);
        }
    }
    Ok(out)
}

/// RFC 7946 lon/lat polygons for a multipolygon (split at the antimeridian).
#[must_use]
pub fn to_lonlat_polygons(mp: &MultiPolygon, opts: &WriteOptions) -> Vec<LonLatPolygon> {
    mp.polygons
        .iter()
        .flat_map(|p| to_rfc7946(p, opts.densify_deg))
        .collect()
}

/// GeoJSON geometry value (`Polygon` when there is exactly one part,
/// `MultiPolygon` otherwise).
#[must_use]
pub fn to_geometry(mp: &MultiPolygon, opts: &WriteOptions) -> Value {
    let parts = to_lonlat_polygons(mp, opts);
    if parts.len() == 1 {
        json!({ "type": "Polygon", "coordinates": parts[0] })
    } else {
        json!({ "type": "MultiPolygon", "coordinates": parts })
    }
}

/// GeoJSON `Feature` with the given properties.
#[must_use]
pub fn to_feature(mp: &MultiPolygon, properties: Map<String, Value>, opts: &WriteOptions) -> Value {
    json!({ "type": "Feature", "geometry": to_geometry(mp, opts), "properties": properties })
}

/// Serialises a multipolygon as a GeoJSON geometry string.
#[must_use]
pub fn to_string(mp: &MultiPolygon, opts: &WriteOptions) -> String {
    to_geometry(mp, opts).to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    const GC: ReadOptions = ReadOptions {
        edges: EdgeInterpretation::GreatCircle,
        dissolve: false,
    };

    #[test]
    fn parse_all_container_types() {
        let poly = r#"{"type":"Polygon","coordinates":[[[0,0],[1,0],[1,1],[0,1],[0,0]]]}"#;
        assert_eq!(parse_str(poly, &GC).unwrap().polygons.len(), 1);
        let multi =
            r#"{"type":"MultiPolygon","coordinates":[[[[0,0],[1,0],[1,1],[0,0]]],[[[5,5],[6,5],[6,6],[5,5]]]]}"#;
        assert_eq!(parse_str(multi, &GC).unwrap().polygons.len(), 2);
        let fc = format!(
            r#"{{"type":"FeatureCollection","features":[{{"type":"Feature","properties":{{}},"geometry":{poly}}},{{"type":"Feature","geometry":null}}]}}"#
        );
        assert_eq!(parse_str(&fc, &GC).unwrap().polygons.len(), 1);
        let gc = format!(r#"{{"type":"GeometryCollection","geometries":[{poly},{multi}]}}"#);
        assert_eq!(parse_str(&gc, &GC).unwrap().polygons.len(), 3);
        // Altitude is ignored.
        let alt = r#"{"type":"Polygon","coordinates":[[[0,0,10],[1,0,10],[1,1,10],[0,0,10]]]}"#;
        assert_eq!(parse_str(alt, &GC).unwrap().polygons.len(), 1);
    }

    #[test]
    fn parse_errors() {
        let cases = [
            "not json",
            r#"{"coordinates":[]}"#,
            r#"{"type":"Point","coordinates":[0,0]}"#,
            r#"{"type":"Polygon"}"#,
            r#"{"type":"Polygon","coordinates":5}"#,
            r#"{"type":"Polygon","coordinates":[]}"#,
            r#"{"type":"Polygon","coordinates":[5]}"#,
            r#"{"type":"Polygon","coordinates":[[5]]}"#,
            r#"{"type":"Polygon","coordinates":[[[0]]]}"#,
            r#"{"type":"Polygon","coordinates":[[["a",0]]]}"#,
            r#"{"type":"Polygon","coordinates":[[[0,"b"]]]}"#,
            r#"{"type":"MultiPolygon","coordinates":5}"#,
            r#"{"type":"FeatureCollection"}"#,
            r#"{"type":"GeometryCollection"}"#,
        ];
        for c in cases {
            assert!(matches!(parse_str(c, &GC), Err(GeoError::GeoJson(_))), "{c}");
        }
        let degenerate = r#"{"type":"Polygon","coordinates":[[[0,0],[1,0],[0,0]]]}"#;
        assert!(matches!(
            parse_str(degenerate, &GC),
            Err(GeoError::DegeneratePolygon(_))
        ));
        assert!(matches!(
            densify_lonlat(&[[0.0, 0.0]], 0.0),
            Err(GeoError::InvalidArgument(_))
        ));
    }

    #[test]
    fn rfc_polar_cap_input() {
        // RFC 7946 style polar cap: full parallel at 80°N closed via the pole.
        let s = r#"{"type":"Polygon","coordinates":[[[-180,80],[180,80],[180,90],[-180,90],[-180,80]]]}"#;
        let mp = parse_str(s, &ReadOptions::default()).unwrap();
        assert!(mp.contains_lonlat([0.0, 85.0]));
        assert!(mp.contains_lonlat([179.0, 81.0]));
        assert!(!mp.contains_lonlat([0.0, 75.0]));
        // Lon/lat-linear edges: area of the cap above 80°.
        let cap = std::f64::consts::TAU * (1.0 - 10f64.to_radians().cos());
        assert_abs_diff_eq!(mp.area_sr() / cap, 1.0, epsilon = 1e-3);
    }

    #[test]
    fn lonlat_linear_vs_great_circle() {
        let s = r#"{"type":"Polygon","coordinates":[[[0,40],[60,40],[60,50],[0,50],[0,40]]]}"#;
        let lin = parse_str(s, &ReadOptions::default()).unwrap();
        let gc = parse_str(s, &GC).unwrap();
        // The 40°N parallel is south of the great circle between its ends.
        assert!(lin.contains_lonlat([30.0, 41.0]));
        assert!(!gc.contains_lonlat([30.0, 41.0]));
    }

    #[test]
    fn antimeridian_round_trip() {
        let s = r#"{"type":"Polygon","coordinates":[[[170,-10],[-170,-10],[-170,10],[170,10],[170,-10]]]}"#;
        let mp = parse_str(s, &GC).unwrap();
        assert!(mp.contains_lonlat([180.0, 0.0]));
        let out = to_geometry(&mp, &WriteOptions::default());
        assert_eq!(out["type"], "MultiPolygon");
        let back = parse_value(
            &out,
            &ReadOptions {
                edges: EdgeInterpretation::GreatCircle,
                dissolve: true,
            },
        )
        .unwrap();
        assert_eq!(back.polygons.len(), 1);
        assert_abs_diff_eq!(back.area_sr(), mp.area_sr(), epsilon = 1e-9);
        // Split RFC input dissolves into one polygon.
        let split = r#"{"type":"MultiPolygon","coordinates":[[[[170,-10],[180,-10],[180,10],[170,10],[170,-10]]],[[[-180,-10],[-170,-10],[-170,10],[-180,10],[-180,-10]]]]}"#;
        let d = parse_str(
            split,
            &ReadOptions {
                edges: EdgeInterpretation::default(),
                dissolve: true,
            },
        )
        .unwrap();
        assert_eq!(d.polygons.len(), 1);
    }

    #[test]
    fn writers() {
        let mp = parse_str(
            r#"{"type":"Polygon","coordinates":[[[0,0],[1,0],[1,1],[0,1],[0,0]]]}"#,
            &GC,
        )
        .unwrap();
        let s = to_string(&mp, &WriteOptions::default());
        assert!(s.contains("\"Polygon\""));
        let mut props = Map::new();
        props.insert("name".into(), json!("beam-1"));
        let f = to_feature(&mp, props, &WriteOptions { densify_deg: Some(0.1) });
        assert_eq!(f["properties"]["name"], "beam-1");
        assert!(f["geometry"]["coordinates"][0].as_array().unwrap().len() > 30);
        let empty = to_geometry(&MultiPolygon::default(), &WriteOptions::default());
        assert_eq!(empty["type"], "MultiPolygon");
        assert_eq!(to_lonlat_polygons(&mp, &WriteOptions::default()).len(), 1);
    }

    #[test]
    fn densify_reads_exact_half_turn_as_eastward() {
        let out = densify_lonlat(&[[0.0, 0.0], [180.0, 0.0], [180.0, 10.0], [0.0, 10.0]], 10.0).unwrap();
        assert!(out.contains(&[90.0, 0.0]), "edge 0→180 must go east through 90");
        assert!(!out.iter().any(|p| p[0] < 0.0));
    }

    #[test]
    fn densify_reads_long_jump_as_antimeridian_crossing() {
        let out = densify_lonlat(&[[170.0, 0.0], [-170.0, 0.0], [-170.0, 1.0], [170.0, 1.0]], 5.0).unwrap();
        assert!(out.contains(&[175.0, 0.0]));
        assert!(
            !out.iter().any(|p| p[0] > 0.0 && p[0] < 170.0),
            "must not go west through 0"
        );
    }

    #[test]
    fn densify_rejects_nan_step() {
        assert!(matches!(
            densify_lonlat(&[[0.0, 0.0]], f64::NAN),
            Err(GeoError::InvalidArgument(_))
        ));
    }

    #[test]
    fn errors_inside_collections_are_reported() {
        let cases = [
            r#"{"type":"MultiPolygon"}"#,
            r#"{"type":"MultiPolygon","coordinates":[[[[0,0],[1,0],[0,0]]]]}"#,
            r#"{"type":"Feature","geometry":{"type":"Point","coordinates":[0,0]}}"#,
            r#"{"type":"FeatureCollection","features":[{"type":"Feature","geometry":{"type":"LineString","coordinates":[]}}]}"#,
            r#"{"type":"GeometryCollection","geometries":[{"type":"Point","coordinates":[0,0]}]}"#,
        ];
        for c in cases {
            assert!(parse_str(c, &GC).is_err(), "{c}");
        }
        assert!(matches!(parse_str(cases[2], &GC), Err(GeoError::GeoJson(_))));
        assert!(matches!(parse_str(cases[1], &GC), Err(GeoError::DegeneratePolygon(_))));
    }

    #[test]
    fn invalid_densification_step_is_rejected() {
        let s = r#"{"type":"Polygon","coordinates":[[[0,0],[1,0],[1,1],[0,0]]]}"#;
        let opts = ReadOptions {
            edges: EdgeInterpretation::LonLatLinear { max_step_deg: 0.0 },
            dissolve: false,
        };
        assert!(matches!(parse_str(s, &opts), Err(GeoError::InvalidArgument(_))));
    }
}
