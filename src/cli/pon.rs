//! CLI subcommand for Panel-of-Normals (PoN) operations.
//!
//! Provides two sub-subcommands:
//!
//! - `kmerdet pon build` — Build a PoN from normal sample detection results.
//! - `kmerdet pon filter` — Filter detection results against a PoN.

use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

use anyhow::{Context, Result};

use crate::filter::pon::PanelOfNormals;

#[derive(clap::Args, Debug)]
pub struct PonArgs {
    #[command(subcommand)]
    pub command: PonCommand,
}

#[derive(clap::Subcommand, Debug)]
pub enum PonCommand {
    /// Build a panel-of-normals from normal sample detection results
    Build(PonBuildArgs),
    /// Filter detection results against a panel-of-normals
    Filter(PonFilterArgs),
}

#[derive(clap::Args, Debug)]
pub struct PonBuildArgs {
    /// Detection result TSV files from normal samples
    #[arg(long, num_args = 1..)]
    pub normals: Vec<PathBuf>,

    /// Output PoN file (JSON)
    #[arg(short, long)]
    pub output: PathBuf,

    /// Minimum VAF to include a variant in the PoN (default: 0.0)
    #[arg(long, default_value = "0.0")]
    pub min_vaf: f64,
}

#[derive(clap::Args, Debug)]
pub struct PonFilterArgs {
    /// Detection result TSV to filter
    #[arg(long)]
    pub input: PathBuf,

    /// Panel-of-normals file (JSON)
    #[arg(long)]
    pub pon: PathBuf,

    /// Frequency threshold: filter variants seen in >= this fraction of normals (default: 0.05 = 5%)
    #[arg(long, default_value = "0.05")]
    pub frequency: f64,

    /// Output file (default: stdout)
    #[arg(short, long)]
    pub output: Option<PathBuf>,
}

pub fn run(args: PonArgs, global: &super::GlobalOptions) -> Result<()> {
    match args.command {
        PonCommand::Build(build_args) => run_build(build_args),
        PonCommand::Filter(filter_args) => run_filter(filter_args, global),
    }
}

fn run_build(args: PonBuildArgs) -> Result<()> {
    // Validate inputs.
    anyhow::ensure!(
        !args.normals.is_empty(),
        "at least one normal sample file is required (--normals)"
    );
    for path in &args.normals {
        anyhow::ensure!(
            path.exists(),
            "normal sample file not found: {}",
            path.display()
        );
    }
    anyhow::ensure!(
        args.min_vaf >= 0.0,
        "--min-vaf must be >= 0.0, got {}",
        args.min_vaf
    );

    let file_refs: Vec<&std::path::Path> = args.normals.iter().map(|p| p.as_path()).collect();

    tracing::info!(
        "building PoN from {} normal sample file(s)",
        file_refs.len()
    );

    let pon = PanelOfNormals::build(&file_refs, args.min_vaf)?;

    tracing::info!(
        "PoN built: {} unique variants from {} samples",
        pon.len(),
        pon.total_samples
    );

    pon.save(&args.output)?;

    tracing::info!("PoN saved to {}", args.output.display());

    // Also print a summary to stderr for the user.
    eprintln!(
        "Panel-of-Normals built: {} variants from {} samples -> {}",
        pon.len(),
        pon.total_samples,
        args.output.display()
    );

    Ok(())
}

fn run_filter(args: PonFilterArgs, global: &super::GlobalOptions) -> Result<()> {
    // Validate inputs.
    anyhow::ensure!(
        args.input.exists(),
        "input file not found: {}",
        args.input.display()
    );
    anyhow::ensure!(
        args.pon.exists(),
        "PoN file not found: {}",
        args.pon.display()
    );
    anyhow::ensure!(
        (0.0..=1.0).contains(&args.frequency),
        "--frequency must be between 0.0 and 1.0, got {}",
        args.frequency
    );

    tracing::info!("loading PoN from {}", args.pon.display());
    let pon = PanelOfNormals::load(&args.pon)?;
    tracing::info!(
        "PoN loaded: {} variants from {} samples (created {})",
        pon.len(),
        pon.total_samples,
        pon.creation_date
    );

    tracing::info!("parsing detection results from {}", args.input.display());
    let calls = crate::output::parse_detection_tsv(&args.input)?;
    let total_calls = calls.len();

    let (passed, filtered) = pon.filter_calls(calls, args.frequency);

    tracing::info!(
        "PoN filter: {}/{} calls passed, {} filtered (frequency threshold: {:.1}%)",
        passed.len(),
        total_calls,
        filtered.len(),
        args.frequency * 100.0
    );

    // Determine output destination.
    let output_path = args.output.as_ref().or(global.output.as_ref());
    let mut output: Box<dyn Write> = match output_path {
        Some(path) => Box::new(
            File::create(path)
                .with_context(|| format!("creating output file: {}", path.display()))?,
        ),
        None => Box::new(std::io::stdout().lock()),
    };

    // Write passed calls in the original detection TSV format.
    crate::output::write_calls(&passed, &global.format, &mut output, global.no_header)?;

    // Print summary to stderr.
    eprintln!(
        "PoN filter: {}/{} calls passed, {} removed (freq >= {:.1}%)",
        passed.len(),
        total_calls,
        filtered.len(),
        args.frequency * 100.0
    );

    // If verbose, list filtered variants.
    if global.verbose >= 1 {
        for call in &filtered {
            if let Some(entry) = pon.check(call, args.frequency) {
                eprintln!(
                    "  PON_FILTERED: {} {}:{} {}/{} (seen in {}/{} normals, freq={:.1}%, max_vaf={:.4})",
                    call.target,
                    call.chrom.as_deref().unwrap_or("?"),
                    call.pos.map(|p| p.to_string()).unwrap_or_else(|| "?".to_string()),
                    call.ref_allele.as_deref().unwrap_or("?"),
                    call.alt_allele.as_deref().unwrap_or("?"),
                    entry.sample_count,
                    entry.total_samples,
                    entry.frequency * 100.0,
                    entry.max_vaf,
                );
            }
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pon_build_args_defaults() {
        // Verify that default values parse correctly via clap.
        use clap::Parser;

        #[derive(Parser)]
        struct TestCli {
            #[command(flatten)]
            args: PonBuildArgs,
        }

        let cli = TestCli::parse_from(["test", "--normals", "a.tsv", "b.tsv", "-o", "pon.json"]);
        assert_eq!(cli.args.normals.len(), 2);
        assert_eq!(cli.args.output, PathBuf::from("pon.json"));
        assert!((cli.args.min_vaf - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_pon_filter_args_defaults() {
        use clap::Parser;

        #[derive(Parser)]
        struct TestCli {
            #[command(flatten)]
            args: PonFilterArgs,
        }

        let cli = TestCli::parse_from([
            "test",
            "--input",
            "results.tsv",
            "--pon",
            "pon.json",
        ]);
        assert_eq!(cli.args.input, PathBuf::from("results.tsv"));
        assert_eq!(cli.args.pon, PathBuf::from("pon.json"));
        assert!((cli.args.frequency - 0.05).abs() < 1e-10);
        assert!(cli.args.output.is_none());
    }

    #[test]
    fn test_pon_filter_args_custom_frequency() {
        use clap::Parser;

        #[derive(Parser)]
        struct TestCli {
            #[command(flatten)]
            args: PonFilterArgs,
        }

        let cli = TestCli::parse_from([
            "test",
            "--input",
            "results.tsv",
            "--pon",
            "pon.json",
            "--frequency",
            "0.10",
        ]);
        assert!((cli.args.frequency - 0.10).abs() < 1e-10);
    }
}
