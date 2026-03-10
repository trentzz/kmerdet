use anyhow::Result;
use clap::Parser;

fn main() -> Result<()> {
    let cli = kmerdet::cli::Cli::parse();
    kmerdet::run(cli)
}
