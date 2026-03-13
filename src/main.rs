use anyhow::Result;
use clap::{CommandFactory, FromArgMatches};

fn main() -> Result<()> {
    // Parse twice: once for the typed struct, once for the raw matches.
    // The raw matches let us distinguish "user explicitly passed --count 2"
    // from "count defaulted to 2", which is needed for config-file merging.
    let matches = kmerdet::cli::Cli::command().get_matches();
    let cli = kmerdet::cli::Cli::from_arg_matches(&matches)
        .map_err(|e: clap::Error| e.exit())
        .unwrap();
    kmerdet::run(cli, matches)
}
