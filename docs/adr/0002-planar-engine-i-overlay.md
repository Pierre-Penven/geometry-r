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
