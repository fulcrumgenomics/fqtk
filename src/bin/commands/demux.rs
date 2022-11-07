use crate::commands::command::Command;
use clap::Parser;

#[derive(Parser, Debug)]
pub(crate) struct Demux {}

impl Command for Demux {
    fn execute(&self) -> anyhow::Result<()> {
        println!("Hello, world!");
        Ok(())
    }
}
