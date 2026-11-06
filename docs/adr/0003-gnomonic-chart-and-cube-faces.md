# 0003. Gnomonic chart with cube-face fallback for spherical booleans

Date: 2026-11-06
Status: Accepted

## Context

Spherical booleans must reuse a planar engine ([ADR 0002](0002-planar-engine-i-overlay.md))
without distorting great-circle edges, and must also handle regions larger than
a hemisphere (for example geostationary coverage unions).

## Options considered

1. Lon/lat plane with antimeridian splitting. Edges stop being great circles,
   and the poles break it.
2. One gnomonic chart. Exact (great circles become lines), but limited to less
   than a hemisphere.
3. A gnomonic chart when the operands fit, and the six cube faces otherwise.

## Decision

Option 3.
- **Chart:** used when every operand is small and their near-minimal
  enclosing cap (Bădoiu–Clarkson) is ≤ 85°.
- **Cube faces:**
  - Each ring is clipped to the face with Sutherland–Hodgman against four
    great-circle planes.
  - The winding offset is fixed per face with one exact point-in-polygon test.
  - Each face is overlaid separately.
- **Merge:** pieces are merged in one chart when the result fits.

## Consequences

- Beams, service areas and countries take the exact, fast chart path.
- Results larger than a hemisphere are returned as cube-face pieces, with
  seams. Area and membership are correct. Seam-free stitching is roadmap work.
- Correctness depends on robust point-in-polygon (half-open arc-crossing
  rule).
