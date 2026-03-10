use std::path::PathBuf;

use anyhow::Result;

#[derive(clap::Args, Debug)]
pub struct DetectArgs {
    /// Jellyfish database file
    #[arg(short = 'd', long)]
    pub db: PathBuf,

    /// Target FASTA files or directories (recursive)
    #[arg(short = 'T', long, num_args = 1..)]
    pub targets: Vec<PathBuf>,

    /// Override k-mer length [default: from DB header]
    #[arg(short = 'k', long)]
    pub kmer_length: Option<u8>,

    /// Min absolute k-mer count
    #[arg(long, default_value = "2")]
    pub count: u32,

    /// Min count ratio threshold
    #[arg(long, default_value = "0.05")]
    pub ratio: f64,

    /// DFS max stack depth
    #[arg(long, default_value = "500")]
    pub max_stack: usize,

    /// DFS max break count
    #[arg(long, default_value = "10")]
    pub max_break: usize,

    /// Max graph nodes
    #[arg(long, default_value = "10000")]
    pub max_node: usize,

    /// Enable overlapping mutation clustering
    #[arg(long)]
    pub cluster: bool,
}

pub fn run(_args: DetectArgs, _global: &super::GlobalOptions) -> Result<()> {
    anyhow::bail!("detect subcommand not yet implemented (Phase 6)")
}
