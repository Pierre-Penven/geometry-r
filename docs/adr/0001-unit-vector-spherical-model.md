# 0001. Spherical model on unit vectors with great-circle edges

Date: 2026-11-06
Status: Accepted

## Context

The library must be correct across the antimeridian and at the poles. Working
in lon/lat makes ±180° and ±90° special cases everywhere: unwrapping, polar
singularities, "which way round" ambiguities.

## Options considered

1. Lon/lat planar geometry with antimeridian special cases (the typical GIS
   approach).
2. Ellipsoidal geodesic polygons.
3. Spherical geometry on 3D unit vectors, with great-circle edges (as in S2).

## Decision

Option 3.
- Points become unit vectors built from geodetic `[lon, lat]`, with latitude
  treated as spherical, as in GeoJSON and S2.
- Edges are minor great-circle arcs.
- Each ring bounds the region on its left.
- Longitude appears only at the API boundaries.

## Consequences

- No antimeridian or pole special cases in the core algorithms.
- Areas are spherical (authalic radius). They differ slightly from
  ellipsoidal areas.
- GeoJSON input with lon/lat-linear edges must be densified on reading. This
  is the default in `geojson::EdgeInterpretation`.
- Ring orientation must be decided explicitly: `Orientation::SmallestArea`
  or `AsGiven`.
