# 0002. i_overlay as the planar boolean engine

Date: 2026-11-06
Status: Accepted

## Context

Robust polygon booleans are hard to write. We need a proven planar engine that
builds on Rust 1.81 ([ADR 0004](0004-pin-rust-1-81.md)).

## Options considered

1. `geo` (BooleanOps). Version 0.29 pulls in `spade`, which needs a newer rustc,
   so it does not build.
2. `i_overlay` 1.10. It builds on 1.81, is fast, and supports fill rules
   including `Positive`.
3. A custom spherical Martinez or Greiner–Hormann. Too costly and
   error-prone.

## Decision

Use `i_overlay` 1.10 directly, through `SingleFloatOverlay::overlay` with
`FillRule::Positive`.

## Consequences

- i_overlay fills **clockwise** contours under `Positive`, and its output
  outer contours are clockwise. We feed mirrored-y coordinates and re-orient
  every output in one adapter (`polygon::planar::orient_shape`).
- Floats are snapped to an adaptive integer grid (about extent/2³⁰).
- Vertices that are collinear in the plane are removed. This is harmless in a
  gnomonic chart, but the RFC 7946 output path must re-insert them
  (`antimeridian::Restorer`).

## Update (2026-11-06): i_overlay 4.5

- Moving to Rust 1.88 ([ADR 0005](0005-edition-2024-msrv-1-88.md)) allowed
  `i_overlay` 4.5, which also removes advisory RUSTSEC-2025-0165 (`i_tree`).
- In 4.x, `FillRule::Positive` fills **counter-clockwise** contours and output
  exteriors are counter-clockwise: the 1.x convention described above was
  reversed. The y-mirroring workaround was removed; `orient_shape` still
  checks every output.
- Collinear-vertex removal is unchanged, so `antimeridian::Restorer` stays.
- `geo` would now build, but i_overlay is kept: it is what `geo` uses
  internally, and depending on it directly keeps the dependency tree small.
