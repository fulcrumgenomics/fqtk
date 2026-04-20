//! `fqtk` — the FASTQ toolkit CLI entry point. Dispatches to subcommand implementations
//! (`demux`, `subsample`, `trim`) via the [`Command`] trait and `enum_dispatch`. Installs
//! `mimalloc` as the global allocator and, on x86_64, runs a startup CPU check against
//! the AVX2 baseline that `.cargo/config.toml` requires.

extern crate core;

pub mod commands;

use anyhow::Result;
use clap::Parser;
use commands::command::Command;
use commands::demux::Demux;
use commands::subsample::Subsample;
use commands::trim::Trim;
use enum_dispatch::enum_dispatch;
use env_logger::Env;

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// Top-level `fqtk` CLI parser. Holds a single [`Subcommand`] variant chosen by the
/// user's argv; every subcommand enum variant in turn wraps that subcommand's own
/// clap-derived options struct.
#[derive(Parser, Debug)]
struct Args {
    #[clap(subcommand)]
    subcommand: Subcommand,
}

/// Exhaustive list of `fqtk` subcommands, wired to the [`Command`] trait via
/// `enum_dispatch` so `execute()` forwards to the chosen variant without a match arm
/// per-subcommand.
#[enum_dispatch(Command)]
#[derive(Parser, Debug)]
#[command(version)]
// Single enum instance per program; the size asymmetry between Trim (large CLI
// struct) and Demux/Subsample is not worth a Box layer on the hot path.
#[allow(clippy::large_enum_variant)]
enum Subcommand {
    Demux(Demux),
    Subsample(Subsample),
    Trim(Trim),
}

/// Best-effort guard against running an AVX2-compiled binary on a CPU that lacks
/// AVX2 (pre-2013 Intel Haswell / pre-2017 AMD Zen / pre-2021 Intel Gracemont
/// Atom). The crate's `.cargo/config.toml` targets `x86-64-v3` for x86_64
/// release builds, which emits AVX2 instructions throughout the binary; hitting
/// one of those on an older CPU produces `SIGILL` with no explanation.
///
/// This runs at the top of `main()` and probes CPUID directly (the `cpuid`
/// instruction carries no SIMD baggage, so the check itself is safe on any
/// x86_64 CPU). If AVX2 is missing we print a clear message and exit 1.
///
/// **Caveat:** the guard only covers code that runs after `main()` starts. If
/// the Rust runtime emits AVX2 ops during pre-main startup (TLS setup, allocator
/// init, etc.), the `SIGILL` will still beat us. In practice Rust's pre-main
/// path is small enough that we don't observe this, but we can't guarantee it.
#[cfg(target_arch = "x86_64")]
fn ensure_avx2_or_die() {
    // `__cpuid`/`__cpuid_count` were stabilized as safe in Rust 1.89 — the `cpuid`
    // instruction is unconditionally available on every x86_64 CPU and has no memory
    // side effects.
    use std::arch::x86_64::{__cpuid, __cpuid_count};
    // Leaf 1: ECX bit 27 = OSXSAVE (OS supports xsave, required to use YMM state);
    // bit 28 = AVX. Leaf 7 sub-leaf 0: EBX bit 5 = AVX2.
    let l1 = __cpuid(1);
    let osxsave = (l1.ecx >> 27) & 1 != 0;
    let avx = (l1.ecx >> 28) & 1 != 0;
    let l7 = __cpuid_count(7, 0);
    let avx2 = (l7.ebx >> 5) & 1 != 0;
    if !(osxsave && avx && avx2) {
        eprintln!(
            "error: this fqtk binary was built for x86-64-v3 (AVX2+) but this\n\
             CPU does not report AVX2 support. Required features: AVX, AVX2,\n\
             OSXSAVE. Rebuild from source with a portable baseline:\n\
             \n\
             \tRUSTFLAGS=\"-C target-cpu=x86-64\" cargo build --release\n"
        );
        std::process::exit(1);
    }
}

/// Process entry point. Runs the AVX2 guard (x86_64 only), initializes env_logger with
/// a default of `info` level, parses argv, and dispatches to the selected subcommand's
/// `execute()`. Errors propagate out as `anyhow::Error` for the runtime's default
/// `eprintln!` + nonzero-exit handling.
fn main() -> Result<()> {
    #[cfg(target_arch = "x86_64")]
    ensure_avx2_or_die();
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();
    let args: Args = Args::parse();
    args.subcommand.execute()
}
