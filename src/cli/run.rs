use std::path::PathBuf;

use anyhow::{Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;

use crate::jellyfish::KmerDatabase;
use crate::sequence::path::KmerPath;
use crate::sequence::target::{RefSeq, Target};
use crate::variant::classifier::Classification;
use crate::variant::quantifier::Quantification;
use crate::variant::{VariantCall, VariantType};
use crate::walker::WalkerConfig;

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

/// Open a jellyfish database (only available when compiled with jellyfish support).
#[cfg(has_jellyfish)]
fn open_db(path: &std::path::Path) -> Result<Box<dyn KmerDatabase>> {
    let db = crate::jellyfish::ffi::JellyfishDb::open(path)?;
    Ok(Box::new(db))
}

/// Stub when compiled without jellyfish support.
#[cfg(not(has_jellyfish))]
fn open_db(_path: &std::path::Path) -> Result<Box<dyn KmerDatabase>> {
    anyhow::bail!(
        "kmerdet was compiled without jellyfish support.\n\
         Install jellyfish 2.2+ and rebuild, or provide a pre-counted database."
    )
}

/// Determine the k-mer length from CLI args or the DB header.
fn resolve_kmer_length(db_path: &std::path::Path) -> Result<u8> {
    let header = crate::jellyfish::header::JfHeader::from_file(db_path)
        .context("reading jellyfish header to determine k-mer length")?;

    header
        .kmer_length()
        .context("could not determine k-mer length from DB header; use --kmer-length to specify")
}

/// Derive a sample name from the database file path.
fn sample_name(db_path: &std::path::Path) -> String {
    db_path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown")
        .to_string()
}

/// Process a single target: walk, build graph, find paths, classify, quantify.
fn process_target(
    target: &Target,
    db: &dyn KmerDatabase,
    k: u8,
    walker_config: &WalkerConfig,
    sample: &str,
) -> Result<Vec<VariantCall>> {
    // a. Decompose target into reference k-mers.
    let refseq = RefSeq::from_target(target.clone(), k)?;

    // b. Walk the k-mer graph starting from reference k-mers.
    let walk_result = crate::walker::walk(db, &refseq.kmers, walker_config);

    // c. Build the directed weighted graph.
    let graph = crate::graph::builder::build_graph(&walk_result, &refseq.kmers);

    // d. Find alternative paths through the graph.
    let paths = crate::graph::pathfind::find_alternative_paths(&graph);

    if paths.is_empty() {
        // No paths at all (disconnected graph) -- produce a Reference call anyway.
        return Ok(vec![make_reference_call(sample, target, db, &refseq)]);
    }

    // The first path is always the reference path.
    let ref_path = &paths[0];

    if paths.len() == 1 {
        // Only reference path: no variants detected.
        return Ok(vec![make_reference_call(sample, target, db, &refseq)]);
    }

    // e. Quantify all paths together via NNLS.
    let quant = crate::variant::quantifier::quantify(&paths, db);

    // f. Classify and build VariantCall for each alternative path.
    let mut calls: Vec<VariantCall> = Vec::new();

    // Include the reference call at index 0.
    calls.push(build_variant_call(
        sample,
        &target.name,
        &Classification {
            variant_type: VariantType::Reference,
            variant_name: String::new(),
            ref_allele: String::new(),
            alt_allele: String::new(),
            start: 0,
            end: 0,
        },
        ref_path,
        ref_path,
        &quant,
        0,
    ));

    // Process each alternative path (indices 1..).
    for (i, alt_path) in paths.iter().enumerate().skip(1) {
        let classification = crate::variant::classifier::classify(ref_path, alt_path, k);

        calls.push(build_variant_call(
            sample,
            &target.name,
            &classification,
            ref_path,
            alt_path,
            &quant,
            i,
        ));
    }

    Ok(calls)
}

/// Build a VariantCall from classification and quantification results.
fn build_variant_call(
    sample: &str,
    target_name: &str,
    classification: &Classification,
    ref_path: &KmerPath,
    alt_path: &KmerPath,
    quant: &Quantification,
    path_index: usize,
) -> VariantCall {
    VariantCall {
        sample: sample.to_string(),
        target: target_name.to_string(),
        variant_type: classification.variant_type,
        variant_name: classification.variant_name.clone(),
        rvaf: quant.rvafs[path_index],
        expression: quant.coefficients[path_index],
        min_coverage: quant.min_coverages[path_index],
        start_kmer_count: 0, // TODO: first k-mer count from path
        ref_sequence: ref_path.to_sequence(),
        alt_sequence: alt_path.to_sequence(),
        info: "vs_ref".to_string(),
        chrom: None,
        pos: None,
        ref_allele: Some(classification.ref_allele.clone()),
        alt_allele: Some(classification.alt_allele.clone()),
    }
}

/// Create a Reference VariantCall for a target with no detected variants.
fn make_reference_call(
    sample: &str,
    target: &Target,
    db: &dyn KmerDatabase,
    refseq: &RefSeq,
) -> VariantCall {
    // Compute mean reference k-mer count as the expression value.
    let counts: Vec<u64> = refseq.kmers.iter().map(|kmer| db.query(kmer)).collect();
    let mean_count = if counts.is_empty() {
        0.0
    } else {
        counts.iter().sum::<u64>() as f64 / counts.len() as f64
    };
    let min_count = counts.iter().copied().min().unwrap_or(0);

    VariantCall {
        sample: sample.to_string(),
        target: target.name.clone(),
        variant_type: VariantType::Reference,
        variant_name: String::new(),
        rvaf: 1.0,
        expression: mean_count,
        min_coverage: min_count,
        start_kmer_count: counts.first().copied().unwrap_or(0),
        ref_sequence: target.sequence.clone(),
        alt_sequence: target.sequence.clone(),
        info: "vs_ref".to_string(),
        chrom: None,
        pos: None,
        ref_allele: None,
        alt_allele: None,
    }
}

/// Parse expected variants from a TSV file.
///
/// Expected format: CHROM\tPOS\tREF\tALT\tTYPE (tab-separated, with header)
fn parse_expected_variants(
    path: &std::path::Path,
) -> Result<Vec<crate::filter::ExpectedVariant>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("opening expected variants: {}", path.display()))?;

    let mut variants = Vec::new();
    for record in reader.records() {
        let record = record?;
        variants.push(crate::filter::ExpectedVariant {
            chrom: record.get(0).unwrap_or("").to_string(),
            pos: record
                .get(1)
                .unwrap_or("0")
                .parse()
                .with_context(|| "parsing POS column in expected variants")?,
            ref_allele: record.get(2).unwrap_or("").to_string(),
            alt_allele: record.get(3).unwrap_or("").to_string(),
            variant_type: record.get(4).unwrap_or("").to_string(),
        });
    }
    Ok(variants)
}

/// Column headers for filter output.
const FILTER_HEADERS: &[&str] = &[
    "sample",
    "chrom",
    "pos",
    "ref_allele",
    "alt_allele",
    "variant_type",
    "found",
    "filter_notes",
    "kmer_vaf",
    "kmer_min_coverage",
    "kmer_expression",
    "ref_sequence",
    "variant_sequence",
];

/// Write filter results as TSV.
fn write_filter_results(
    results: &[crate::filter::FilterResult],
    writer: &mut dyn std::io::Write,
    no_header: bool,
) -> Result<()> {
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(writer);

    if !no_header {
        wtr.write_record(FILTER_HEADERS)?;
    }

    for r in results {
        wtr.write_record(&[
            &r.sample,
            &r.chrom,
            &r.pos.to_string(),
            &r.ref_allele,
            &r.alt_allele,
            &r.variant_type,
            &r.found,
            &r.filter_notes,
            &r.kmer_vaf
                .map(|v| format!("{:.6}", v))
                .unwrap_or_default(),
            &r.kmer_min_coverage
                .map(|v| v.to_string())
                .unwrap_or_default(),
            &r.kmer_expression
                .map(|v| format!("{:.2}", v))
                .unwrap_or_default(),
            &r.ref_sequence.clone().unwrap_or_default(),
            &r.variant_sequence.clone().unwrap_or_default(),
        ])?;
    }

    wtr.flush()?;
    Ok(())
}

pub fn run(args: RunArgs, global: &super::GlobalOptions) -> Result<()> {
    use std::time::Instant;

    let total_start = Instant::now();

    // 1. Open jellyfish database.
    tracing::info!("Opening jellyfish database: {}", args.db.display());
    let db = open_db(&args.db)?;
    let db_name = sample_name(&args.db);

    // 2. Determine k-mer length from DB header.
    let k = resolve_kmer_length(&args.db)?;
    tracing::info!("Using k-mer length: {}", k);

    // 3. Load all target FASTA files from the targets directory.
    let targets = crate::sequence::target::load_targets(&args.targets_dir)
        .with_context(|| format!("loading targets from {}", args.targets_dir.display()))?;
    tracing::info!("Loaded {} targets", targets.len());

    if targets.is_empty() {
        anyhow::bail!("no targets loaded from {}; check --targets-dir path", args.targets_dir.display());
    }

    // ── Stage 1: Detection ──────────────────────────────────────────────
    tracing::info!("Stage 1/2: Detection...");
    let stage1_start = Instant::now();

    let walker_config = WalkerConfig {
        count: args.count,
        ratio: args.ratio,
        ..Default::default()
    };

    // Set up progress bar.
    let pb = ProgressBar::new(targets.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} targets ({eta})",
            )
            .unwrap()
            .progress_chars("=>-"),
    );

    // Process each target in parallel.
    let all_calls: Vec<VariantCall> = targets
        .par_iter()
        .flat_map(|target| {
            let result = process_target(target, db.as_ref(), k, &walker_config, &db_name);
            pb.inc(1);
            match result {
                Ok(calls) => calls,
                Err(e) => {
                    tracing::warn!("target '{}' failed: {}", target.name, e);
                    vec![]
                }
            }
        })
        .collect();

    pb.finish_with_message("detection done");

    let stage1_elapsed = stage1_start.elapsed();
    tracing::info!(
        "Detection complete: {} variant calls across {} targets ({:.1}s)",
        all_calls.len(),
        targets.len(),
        stage1_elapsed.as_secs_f64()
    );

    // ── Stage 2: Filtering ──────────────────────────────────────────────
    tracing::info!("Stage 2/2: Filtering...");
    let stage2_start = Instant::now();

    // Parse expected variants from the --expected file.
    let expected = parse_expected_variants(&args.expected)?;
    tracing::info!("Loaded {} expected variants", expected.len());

    // Build FilterConfig from CLI args.
    let filter_config = crate::filter::FilterConfig {
        min_coverage: args.min_coverage,
        min_vaf: args.min_vaf,
        min_expression: 0.0,
        use_alt: args.use_alt,
        types: vec![],
    };

    // Run filtering.
    let results = crate::filter::filter_results(&all_calls, &expected, &filter_config)?;

    let stage2_elapsed = stage2_start.elapsed();

    let found_count = results.iter().filter(|r| r.found == "Found").count();
    tracing::info!(
        "Filtering complete: {}/{} expected variants found ({:.1}s)",
        found_count,
        results.len(),
        stage2_elapsed.as_secs_f64()
    );

    // ── Write output ────────────────────────────────────────────────────
    let mut output: Box<dyn std::io::Write> = match &global.output {
        Some(path) => Box::new(
            std::fs::File::create(path)
                .with_context(|| format!("creating output file: {}", path.display()))?,
        ),
        None => Box::new(std::io::stdout().lock()),
    };

    write_filter_results(&results, &mut output, global.no_header)?;

    let total_elapsed = total_start.elapsed();
    tracing::info!(
        "Pipeline complete ({:.1}s total: detect {:.1}s + filter {:.1}s)",
        total_elapsed.as_secs_f64(),
        stage1_elapsed.as_secs_f64(),
        stage2_elapsed.as_secs_f64()
    );

    Ok(())
}
