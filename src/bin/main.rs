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

fn main() -> Result<()> {
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();
    let args: Args = Args::parse();
    args.subcommand.execute()
}
