use std::path::PathBuf;

use anyhow::{Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;

use crate::graph::pruning::PruneConfig;
use crate::jellyfish::KmerDatabase;
use crate::sequence::path::KmerPath;
use crate::sequence::target::{RefSeq, Target};
use crate::variant::classifier::Classification;
use crate::variant::quantifier::Quantification;
use crate::variant::{VariantCall, VariantType};
use crate::walker::WalkerConfig;

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

    /// Disable graph pruning before pathfinding
    #[arg(long)]
    pub no_prune: bool,

    /// Minimum node count for pruning (0 = adaptive)
    #[arg(long, default_value = "0")]
    pub prune_min_count: u64,

    /// Error rate for adaptive pruning threshold
    #[arg(long, default_value = "0.001")]
    pub prune_error_rate: f64,

    /// Bubble coverage ratio for pruning
    #[arg(long, default_value = "0.10")]
    pub prune_bubble_ratio: f64,
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
fn resolve_kmer_length(args_k: Option<u8>, db_path: &std::path::Path) -> Result<u8> {
    if let Some(k) = args_k {
        return Ok(k);
    }

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

/// Process a single target: walk, build graph, prune, find paths, classify, quantify.
fn process_target(
    target: &Target,
    db: &dyn KmerDatabase,
    k: u8,
    walker_config: &WalkerConfig,
    prune_config: Option<&PruneConfig>,
    sample: &str,
) -> Result<Vec<VariantCall>> {
    // a. Decompose target into reference k-mers.
    let refseq = RefSeq::from_target(target.clone(), k)?;

    // b. Walk the k-mer graph starting from reference k-mers.
    let walk_result = crate::walker::walk(db, &refseq.kmers, walker_config);

    // c. Build the directed weighted graph.
    let mut graph = crate::graph::builder::build_graph(&walk_result, &refseq.kmers);

    // c2. Prune the graph if enabled.
    if let Some(prune_cfg) = prune_config {
        let removed = crate::graph::pruning::prune_graph(&mut graph, prune_cfg, k);
        if removed > 0 {
            tracing::debug!("pruned {} nodes from graph for target '{}'", removed, target.name);
        }
    }

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
        let classification =
            crate::variant::classifier::classify(ref_path, alt_path, k);

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

pub fn run(args: DetectArgs, global: &super::GlobalOptions) -> Result<()> {
    // 1. Open jellyfish database.
    let db = open_db(&args.db)?;
    let db_name = sample_name(&args.db);

    // 2. Determine k-mer length (CLI override > DB header).
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

    // 4. Build walker configuration from args.
    let walker_config = WalkerConfig {
        count: args.count,
        ratio: args.ratio,
        max_stack: args.max_stack,
        max_break: args.max_break,
        max_node: args.max_node,
    };

    // 4b. Build prune configuration (unless --no-prune).
    let prune_config = if args.no_prune {
        None
    } else {
        Some(PruneConfig {
            min_count: args.prune_min_count,
            error_rate: args.prune_error_rate,
            clip_tips: true,
            max_tip_fraction: 0.5,
            collapse_bubbles: true,
            bubble_coverage_ratio: args.prune_bubble_ratio,
        })
    };

    // 5. Set up progress bar.
    let pb = ProgressBar::new(targets.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} targets ({eta})",
            )
            .unwrap()
            .progress_chars("=>-"),
    );

    // 6. Process each target in parallel.
    let all_calls: Vec<VariantCall> = targets
        .par_iter()
        .flat_map(|target| {
            let result = process_target(target, db.as_ref(), k, &walker_config, prune_config.as_ref(), &db_name);
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

    pb.finish_with_message("done");

    tracing::info!(
        "detected {} variant calls across {} targets",
        all_calls.len(),
        targets.len()
    );

    // 7. Write output.
    match &global.output {
        Some(path) => {
            if global.format == super::OutputFormat::Xlsx {
                crate::output::write_calls_xlsx(&all_calls, path)?;
                return Ok(());
            }
            let mut file = std::fs::File::create(path)
                .with_context(|| format!("creating output file: {}", path.display()))?;
            crate::output::write_calls(&all_calls, &global.format, &mut file, global.no_header)?;
        }
        None => {
            let mut stdout = std::io::stdout().lock();
            crate::output::write_calls(
                &all_calls,
                &global.format,
                &mut stdout,
                global.no_header,
            )?;
        }
    }

    Ok(())
}
