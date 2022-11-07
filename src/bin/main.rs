pub mod commands;

use anyhow::Result;
use clap::Parser;
use commands::command::Command;
use commands::demux::Demux;
use enum_dispatch::enum_dispatch;
use env_logger::Env;

#[derive(Parser, Debug)]
struct Args {
    #[clap(subcommand)]
    subcommand: Subcommand,
}

#[enum_dispatch(Command)]
#[derive(Parser, Debug)]
enum Subcommand {
    Demux(Demux),
}

fn main() -> Result<()> {
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();
    let args: Args = Args::parse();
    args.subcommand.execute()
}
