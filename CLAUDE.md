# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

fqtk is a FASTQ toolkit written in Rust. Currently contains a single tool, `fqtk demux`, for demultiplexing FASTQs based on sample barcodes with support for IUPAC ambiguity codes, configurable mismatch tolerances, and multi-threaded compressed output.

## Build & Test Commands

```bash
# Full verification (format + clippy + tests) — run before pushing
bash ci/check.sh

# Individual steps
cargo fmt --all
cargo clippy --all-features --all-targets -- -D warnings
cargo test

# Run a single test
cargo test <test_name>

# Build release binary
cargo build --release
```

CI also runs `src/scripts/precommit.sh` (same checks but with `--locked`). After changing CLI arguments or docstrings, update the README usage block:

```bash
cargo build && bash .github/scripts/update-docs.sh
```

CI will fail if README usage is out of date.

## Toolchain

Pinned to Rust 1.85 via `rust-toolchain.toml`. Format settings: `max_width = 100`, `use_small_heuristics = "max"` (see `rustfmt.toml`).

## Architecture

The crate produces both a library (`fqtk_lib`) and a binary (`fqtk`).

### Library (`src/lib/`)

- `mod.rs` — DNA/IUPAC base encoding/decoding using bitmask representation. Constants `DNA_MASKS` and `IUPAC_MASKS` map bases to bit patterns; `encode()`/`decode()` convert between byte sequences and `BitEnc`.
- `bitenc.rs` — `BitEnc`: a fixed-width bit-packed vector (adapted from rust-bio). Stores bases compactly in `u32` blocks. The `hamming()` method supports asymmetric IUPAC fuzzy matching with early termination.
- `barcode_matching.rs` — `BarcodeMatcher`: matches observed barcode sequences against expected sample barcodes. Uses `BitEnc` hamming distance, applies max-mismatch and min-delta filters, and optionally caches results (ahash `HashMap`).
- `samples.rs` — `Sample` and `SampleGroup`: parse and validate sample metadata from TSV files (via `fgoxide::DelimFile`). Enforces unique sample IDs, unique barcodes, and uniform barcode length.

### Binary (`src/bin/`)

- `main.rs` — CLI entry point using clap derive. Uses `enum_dispatch` for subcommand dispatch and `mimalloc` as the global allocator.
- `commands/command.rs` — `Command` trait (with `enum_dispatch`) that each subcommand implements.
- `commands/demux.rs` — The `Demux` subcommand. Core types:
  - `ReadSet` — a single FASTQ record split into typed segments (template, sample barcode, molecular barcode, cellular barcode) per the read structures.
  - `ReadSetIterator` — reads from multiple FASTQ files in parallel, applying read structures to extract segments.
  - `SampleWriters` — manages per-sample, per-segment-type output writers.
  - `DemuxMetric` — per-sample demux statistics written to `demux-metrics.txt`.
  - Uses `pooled-writer` for multi-threaded BGZF-compressed output and `read-structure` crate for read structure parsing.

### Adding a New Subcommand

1. Create `src/bin/commands/<name>.rs` with a struct deriving `Parser` that implements `Command`.
2. Add the module to `src/bin/commands/mod.rs`.
3. Add the variant to the `Subcommand` enum in `src/bin/main.rs`.

## Code Organization

### Module Ordering

Command modules (`src/bin/commands/*.rs`) must follow this ordering convention:

1. **`use` statements**
2. **Constants and type aliases**
3. **Structs/enums, each immediately followed by all its impl blocks:**
   - CLI options / command struct (the `clap::Parser`-derived struct) + `impl` + `impl Command`
   - Per-run config structs (e.g. `PipelineConfig`) + `impl`
   - Metric / output-row structs + `impl`
   - Helper structs/enums — ordered higher-to-lower level (if A uses B, A comes first), then by importance to the overall implementation
4. **Module-level functions** (functions operating on primitives or external types — SIMD kernels, byte-slice helpers, parsers that don't belong to a type we own)
5. **`#[cfg(test)] mod tests`**

### Impl Block Rules

- Consolidate all inherent methods into one `impl` block per struct; keep each trait impl as a separate block.
- Trait impls go immediately after the struct's own `impl` block.
- Within an impl block, order methods callers-before-callees (higher in the call stack first); constructors (`new`, `from_str`) always come first.
- Functions that naturally belong to a type we own should be methods, not standalone functions.

## Key Dependencies

- `read-structure` — parses read structure strings (e.g. `8B92T`) into typed segments.
- `fgoxide` — Fulcrum Genomics utilities for delimited file I/O.
- `pooled-writer` — thread-pool-based BGZF writer for parallel compression.
- `flate2` with `zlib-ng` feature — faster gzip via C backend (requires C compiler).
- `seq_io` — FASTQ parsing.
- `ahash` — fast hashing for the barcode match cache.
