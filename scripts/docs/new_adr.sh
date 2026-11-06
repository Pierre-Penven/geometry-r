#!/usr/bin/env bash
# Create a numbered architecture decision record in docs/adr/.
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "usage: $0 \"decision title\"" >&2
  exit 2
fi

dir="docs/adr"
mkdir -p "$dir"

next=$(printf "%04d" $(( $(ls "$dir" 2>/dev/null | grep -c '^[0-9]\{4\}-') + 1 )))
slug=$(echo "$1" | tr '[:upper:]' '[:lower:]' | sed -E 's/[^a-z0-9]+/-/g; s/^-|-$//g')
file="$dir/$next-$slug.md"

cat > "$file" <<ADR
# $next. $1

Date: $(date +%F)
Status: Proposed

## Context

## Options considered

## Decision

## Consequences
ADR

echo "Created $file"
