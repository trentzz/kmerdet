use std::path::PathBuf;

use anyhow::Result;

#[derive(clap::Args, Debug)]
pub struct MergeArgs {
    /// Input result files to merge
    #[arg(short, long, num_args = 1..)]
    pub input: Vec<PathBuf>,

    /// Deduplicate rows with identical key columns
    #[arg(long)]
    pub deduplicate: bool,
}

pub fn run(_args: MergeArgs, _global: &super::GlobalOptions) -> Result<()> {
    anyhow::bail!("merge subcommand not yet implemented (Phase 7)")
}
