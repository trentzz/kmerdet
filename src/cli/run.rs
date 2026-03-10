use std::path::PathBuf;

use anyhow::Result;

#[derive(clap::Args, Debug)]
pub struct RunArgs {
    /// Jellyfish database file
    #[arg(short = 'd', long)]
    pub db: PathBuf,

    /// Directory containing target FASTA files
    #[arg(long)]
    pub targets_dir: PathBuf,

    /// Expected variants file for filtering
    #[arg(long)]
    pub expected: PathBuf,

    // -- detect params --
    /// Min absolute k-mer count
    #[arg(long, default_value = "2")]
    pub count: u32,

    /// Min count ratio threshold
    #[arg(long, default_value = "0.05")]
    pub ratio: f64,

    /// Enable overlapping mutation clustering
    #[arg(long)]
    pub cluster: bool,

    // -- filter params --
    /// Match by full ALT_SEQUENCE instead of POS/REF/ALT
    #[arg(long)]
    pub use_alt: bool,

    /// Min k-mer coverage for filtering
    #[arg(long, default_value = "3")]
    pub min_coverage: u32,

    /// Min VAF for filtering
    #[arg(long, default_value = "0.0")]
    pub min_vaf: f64,

    // -- optional pipeline stages --
    /// Also run stats after filtering
    #[arg(long)]
    pub stats: bool,

    /// Also generate plots after filtering
    #[arg(long)]
    pub plot: bool,
}

pub fn run(_args: RunArgs, _global: &super::GlobalOptions) -> Result<()> {
    anyhow::bail!("run subcommand not yet implemented (Phase 8)")
}
