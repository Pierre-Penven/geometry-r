.PHONY: all fmt lint test cov bench docs docs-check check

all: check

fmt:
	cargo fmt --all

lint:
	cargo clippy --all-targets --all-features -- -D warnings

test:
	cargo test --all-features

# Requires: rustup component add llvm-tools-preview
#           cargo install cargo-llvm-cov --version 0.6.21 --locked
cov:
	cargo llvm-cov --all-features --summary-only --show-missing-lines

bench:
	cargo run --release --example bench

docs-check:
	@bash scripts/docs/check_docs.sh
	@python3 scripts/docs/check_links.py .

docs: docs-check
	cargo doc --no-deps --all-features
	@echo "Documentation checks passed"

check: lint test docs
	cargo fmt --all --check
