use std::path::PathBuf;

use anyhow::Result;

#[derive(clap::Args, Debug)]
pub struct StatsArgs {
    /// Input result file(s) to summarize
    #[arg(short, long, num_args = 1..)]
    pub input: Vec<PathBuf>,

    /// Group statistics by this column
    #[arg(long)]
    pub group_by: Option<String>,
}

pub fn run(_args: StatsArgs, _global: &super::GlobalOptions) -> Result<()> {
    anyhow::bail!("stats subcommand not yet implemented (Phase 7)")
}
