use std::io::Write;
use std::path::PathBuf;

use anyhow::{Context, Result};

use crate::jellyfish::KmerDatabase;
use crate::sequence::target::{RefSeq, Target};

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

/// Coverage statistics for a single target.
struct CoverageStats {
    target: String,
    num_kmers: usize,
    mean: f64,
    median: f64,
    min: u64,
    max: u64,
    stdev: f64,
    cv: f64,
}

/// Column headers for coverage output.
const COVERAGE_HEADERS: &[&str] = &[
    "target",
    "num_kmers",
    "mean_coverage",
    "median_coverage",
    "min_coverage",
    "max_coverage",
    "stdev",
    "cv",
];

/// Open a jellyfish database using the pure Rust reader.
fn open_db(path: &std::path::Path) -> Result<Box<dyn KmerDatabase>> {
    let reader = crate::jellyfish::reader::JellyfishReader::open(path)?;
    Ok(Box::new(reader))
}

/// Determine the k-mer length from CLI args or the DB header.
fn resolve_kmer_length(args_k: Option<u8>, db_path: &std::path::Path) -> Result<u8> {
    if let Some(k) = args_k {
        return Ok(k);
    }

    // Use jellyfish-reader to read the header
    let header = jellyfish_reader::FileHeader::read(
        &mut std::fs::File::open(db_path)
            .with_context(|| format!("opening jellyfish DB: {}", db_path.display()))?,
    )
    .context("reading jellyfish header")?;

    header
        .k()
        .map(|k| k as u8)
        .context("could not determine k-mer length from DB header; use --kmer-length to specify")
}

/// Compute coverage statistics for a target against the database.
fn compute_coverage(target: &Target, db: &dyn KmerDatabase, k: u8) -> Result<CoverageStats> {
    let refseq = RefSeq::from_target(target.clone(), k)?;

    let counts: Vec<u64> = refseq.kmers.iter().map(|kmer| db.query(kmer)).collect();
    let num_kmers = counts.len();

    if num_kmers == 0 {
        return Ok(CoverageStats {
            target: target.name.clone(),
            num_kmers: 0,
            mean: 0.0,
            median: 0.0,
            min: 0,
            max: 0,
            stdev: 0.0,
            cv: 0.0,
        });
    }

    let min = counts.iter().copied().min().unwrap_or(0);
    let max = counts.iter().copied().max().unwrap_or(0);

    let sum: u64 = counts.iter().sum();
    let mean = sum as f64 / num_kmers as f64;

    // Standard deviation (population).
    let variance = counts
        .iter()
        .map(|&c| {
            let diff = c as f64 - mean;
            diff * diff
        })
        .sum::<f64>()
        / num_kmers as f64;
    let stdev = variance.sqrt();

    // Coefficient of variation.
    let cv = if mean > 0.0 { stdev / mean } else { 0.0 };

    // Median.
    let mut sorted = counts;
    sorted.sort_unstable();
    let median = if num_kmers % 2 == 0 {
        (sorted[num_kmers / 2 - 1] as f64 + sorted[num_kmers / 2] as f64) / 2.0
    } else {
        sorted[num_kmers / 2] as f64
    };

    Ok(CoverageStats {
        target: target.name.clone(),
        num_kmers,
        mean,
        median,
        min,
        max,
        stdev,
        cv,
    })
}

/// Write coverage statistics as delimited text (TSV or CSV).
fn write_coverage_stats(
    stats: &[CoverageStats],
    writer: &mut dyn Write,
    delimiter: u8,
    no_header: bool,
) -> Result<()> {
    let delim = delimiter as char;

    if !no_header {
        writeln!(writer, "{}", COVERAGE_HEADERS.join(&delim.to_string()))?;
    }

    for s in stats {
        writeln!(
            writer,
            "{1}{0}{2}{0}{3:.2}{0}{4:.1}{0}{5}{0}{6}{0}{7:.2}{0}{8:.4}",
            delim,
            s.target,
            s.num_kmers,
            s.mean,
            s.median,
            s.min,
            s.max,
            s.stdev,
            s.cv,
        )?;
    }

    Ok(())
}

pub fn run(args: CoverageArgs, global: &super::GlobalOptions) -> Result<()> {
    // 1. Open jellyfish database.
    let db = open_db(&args.db)?;

    // 2. Determine k-mer length.
    let k = resolve_kmer_length(args.kmer_length, &args.db)?;
    tracing::info!("using k-mer length: {}", k);

    // 3. Load all target FASTA files.
    let mut targets: Vec<Target> = Vec::new();
    for path in &args.targets {
        let loaded = crate::sequence::target::load_targets(path)
            .with_context(|| format!("loading targets from {}", path.display()))?;
        targets.extend(loaded);
    }
    tracing::info!("loaded {} targets", targets.len());

    if targets.is_empty() {
        anyhow::bail!("no targets loaded; check --targets paths");
    }

    // 4. Compute coverage stats for each target.
    let mut all_stats: Vec<CoverageStats> = Vec::with_capacity(targets.len());
    for target in &targets {
        match compute_coverage(target, db.as_ref(), k) {
            Ok(stats) => all_stats.push(stats),
            Err(e) => {
                tracing::warn!("target '{}' failed: {}", target.name, e);
            }
        }
    }

    // 5. Determine delimiter from output format.
    let delimiter = match global.format {
        super::OutputFormat::Csv => b',',
        _ => b'\t', // Default to TSV for all non-CSV formats
    };

    // 6. Write output.
    match &global.output {
        Some(path) => {
            let mut file = std::fs::File::create(path)
                .with_context(|| format!("creating output file: {}", path.display()))?;
            write_coverage_stats(&all_stats, &mut file, delimiter, global.no_header)?;
        }
        None => {
            let mut stdout = std::io::stdout().lock();
            write_coverage_stats(&all_stats, &mut stdout, delimiter, global.no_header)?;
        }
    }

    Ok(())
}
