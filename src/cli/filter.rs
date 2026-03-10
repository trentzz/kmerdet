use std::path::PathBuf;

use anyhow::Result;

#[derive(clap::Args, Debug)]
pub struct FilterArgs {
    /// Detection results input file (from detect/merge)
    #[arg(short, long)]
    pub input: PathBuf,

    /// Expected variants file (TSV: CHROM, POS, REF, ALT, TYPE)
    #[arg(long)]
    pub targets: PathBuf,

    /// Match by full ALT_SEQUENCE instead of POS/REF/ALT
    #[arg(long)]
    pub use_alt: bool,

    /// Min k-mer coverage
    #[arg(long, default_value = "3")]
    pub min_coverage: u32,

    /// Min variant allele frequency
    #[arg(long, default_value = "0.0")]
    pub min_vaf: f64,

    /// Min expression value
    #[arg(long, default_value = "0.0")]
    pub min_expression: f64,

    /// Variant types to include [default: all]
    #[arg(long, num_args = 1..)]
    pub types: Vec<String>,
}

pub fn run(_args: FilterArgs, _global: &super::GlobalOptions) -> Result<()> {
    anyhow::bail!("filter subcommand not yet implemented (Phase 7)")
}
