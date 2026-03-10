use std::path::PathBuf;

use clap::{Parser, Subcommand, ValueEnum};

pub mod coverage;
pub mod detect;
pub mod filter;
pub mod merge;
pub mod plot;
pub mod run;
pub mod stats;

#[derive(Parser, Debug)]
#[command(
    name = "kmerdet",
    version,
    about = "K-mer based variant detection for targeted liquid biopsy monitoring",
    long_about = "kmerdet detects variants by walking k-mer graphs built from jellyfish databases.\n\n\
                  It replaces the km + kmtools + kam ecosystem with a single, fast Rust binary."
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,

    #[command(flatten)]
    pub global: GlobalOptions,
}

#[derive(clap::Args, Debug, Clone)]
pub struct GlobalOptions {
    /// Config file (TOML)
    #[arg(short, long, global = true)]
    pub config: Option<PathBuf>,

    /// Thread count (0 = all available cores)
    #[arg(short = 't', long, global = true, default_value_t = 0)]
    pub threads: usize,

    /// Increase verbosity (-v, -vv, -vvv)
    #[arg(short, long, global = true, action = clap::ArgAction::Count)]
    pub verbose: u8,

    /// Suppress non-error output
    #[arg(short, long, global = true)]
    pub quiet: bool,

    /// Output file [default: stdout]
    #[arg(short, long, global = true)]
    pub output: Option<PathBuf>,

    /// Output format
    #[arg(short, long, global = true, default_value = "tsv")]
    pub format: OutputFormat,

    /// Omit column headers
    #[arg(long, global = true)]
    pub no_header: bool,
}

#[derive(ValueEnum, Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    Tsv,
    Csv,
    Vcf,
    Json,
    Jsonl,
    Xlsx,
}

#[derive(Subcommand, Debug)]
pub enum Command {
    /// Run variant detection on targets against a jellyfish database
    Detect(detect::DetectArgs),
    /// Apply tumor-informed filtering to detection results
    Filter(filter::FilterArgs),
    /// Merge results from multiple detection runs
    Merge(merge::MergeArgs),
    /// Generate summary statistics from results
    Stats(stats::StatsArgs),
    /// Generate visualization charts
    Plot(plot::PlotArgs),
    /// Report k-mer coverage statistics for a database
    Coverage(coverage::CoverageArgs),
    /// Full pipeline: detect -> merge -> filter (-> optional stats/plot)
    Run(run::RunArgs),
}
