use std::path::PathBuf;

use anyhow::Result;

#[derive(clap::Args, Debug)]
pub struct CoverageArgs {
    /// Jellyfish database file
    #[arg(short = 'd', long)]
    pub db: PathBuf,

    /// Target FASTA files or directories
    #[arg(short = 'T', long, num_args = 1..)]
    pub targets: Vec<PathBuf>,

    /// Override k-mer length [default: from DB header]
    #[arg(short = 'k', long)]
    pub kmer_length: Option<u8>,
}

pub fn run(_args: CoverageArgs, _global: &super::GlobalOptions) -> Result<()> {
    anyhow::bail!("coverage subcommand not yet implemented (Phase 1)")
}
