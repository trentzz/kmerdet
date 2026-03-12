use std::path::PathBuf;

use anyhow::{Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;

use crate::graph::pruning::PruneConfig;
use crate::jellyfish::KmerDatabase;
use crate::sequence::path::KmerPath;
use crate::sequence::target::{RefSeq, Target};
use crate::variant::bootstrap::BootstrapConfig;
use crate::variant::classifier::Classification;
use crate::variant::quantifier::Quantification;
use crate::variant::{VariantCall, VariantType};
use crate::walker::adaptive::{AdaptiveConfig, estimate_thresholds};
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

    /// Write a .dot graph file per target into this directory for debugging.
    #[arg(long, value_name = "DIR")]
    pub debug_graph: Option<PathBuf>,

    /// Enable bootstrap confidence intervals on rVAF estimates.
    #[arg(long)]
    pub bootstrap: bool,

    /// Number of bootstrap replicates (default 1000, only used with --bootstrap)
    #[arg(long, default_value = "1000")]
    pub bootstrap_replicates: usize,

    /// Bootstrap confidence level (default 0.95 for 95% CI)
    #[arg(long, default_value = "0.95")]
    pub bootstrap_confidence: f64,

    /// Bootstrap random seed for reproducibility
    #[arg(long)]
    pub bootstrap_seed: Option<u64>,

    /// Use adaptive thresholds based on sample coverage and error rate.
    #[arg(long)]
    pub adaptive: bool,

    /// Walk both forward and backward from reference k-mers.
    #[arg(long)]
    pub bidirectional: bool,

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

/// Derive a sample name from the database file path.
fn sample_name(db_path: &std::path::Path) -> String {
    db_path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("unknown")
        .to_string()
}

/// Default assumed read length for autocorrelation correction in Fisher's method.
const DEFAULT_READ_LENGTH: usize = 150;

/// Process a single target: walk, build graph, prune, find paths, classify, quantify,
/// and compute statistical confidence.
fn process_target(
    target: &Target,
    db: &dyn KmerDatabase,
    k: u8,
    walker_config: &WalkerConfig,
    prune_config: Option<&PruneConfig>,
    sample: &str,
    debug_graph_dir: Option<&std::path::Path>,
    bootstrap_config: Option<&BootstrapConfig>,
) -> Result<Vec<VariantCall>> {
    // a. Decompose target into reference k-mers.
    let refseq = RefSeq::from_target(target.clone(), k)?;

    // b. Walk the k-mer graph starting from reference k-mers.
    let walk_result = if walker_config.bidirectional {
        crate::walker::walk_bidirectional(db, &refseq.kmers, walker_config)
    } else {
        crate::walker::walk(db, &refseq.kmers, walker_config)
    };

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

    // d'. Optionally write a debug DOT file for this target.
    if let Some(dir) = debug_graph_dir {
        let safe_name: String = target
            .name
            .chars()
            .map(|c| if c.is_alphanumeric() || c == '-' || c == '_' { c } else { '_' })
            .collect();
        let dot_path = dir.join(format!("{}.dot", safe_name));

        let path_refs: Vec<&crate::sequence::path::KmerPath> = paths.iter().collect();
        let highlight = if paths.is_empty() {
            None
        } else {
            Some(path_refs.as_slice())
        };
        let dot_str = crate::graph::dot::export_dot(&graph, &target.name, highlight);

        if let Err(e) = crate::graph::dot::write_dot_file(&dot_str, &dot_path) {
            tracing::warn!(
                "failed to write debug graph for '{}': {}",
                target.name,
                e
            );
        } else {
            tracing::debug!("wrote debug graph: {}", dot_path.display());
        }
    }

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
    let mut quant = crate::variant::quantifier::quantify(&paths, db);

    // e'. Optionally run bootstrap to compute confidence intervals.
    if let Some(bs_config) = bootstrap_config {
        let cis =
            crate::variant::bootstrap::bootstrap_confidence_intervals(&paths, db, bs_config);
        quant.ci_lower = cis.iter().map(|ci| ci.ci_lower).collect();
        quant.ci_upper = cis.iter().map(|ci| ci.ci_upper).collect();
    }

    // e'. Estimate per-target error rate from reference k-mers for confidence scoring.
    let error_rate = crate::confidence::pvalue::estimate_error_rate(db, &refseq.kmers);

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

        // Compute statistical confidence for this variant.
        let pvalue = crate::confidence::compute_variant_pvalue(
            db,
            alt_path,
            ref_path,
            error_rate,
            DEFAULT_READ_LENGTH,
        );
        let qual = crate::confidence::phred_qual(pvalue);

        let mut call = build_variant_call(
            sample,
            &target.name,
            &classification,
            ref_path,
            alt_path,
            &quant,
            i,
        );
        call.pvalue = Some(pvalue);
        call.qual = Some(qual);

        calls.push(call);
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
    let ci_lower = quant.ci_lower.get(path_index).copied();
    let ci_upper = quant.ci_upper.get(path_index).copied();
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
        pvalue: None,
        qual: None,
        ci_lower,
        ci_upper,
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
        pvalue: None,
        qual: None,
        ci_lower: None,
        ci_upper: None,
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
    let (count, ratio) = if args.adaptive {
        // Collect reference k-mers from the first few targets for estimation.
        let mut sample_ref_kmers: Vec<String> = Vec::new();
        let max_targets_to_sample = 10.min(targets.len());
        for target in &targets[..max_targets_to_sample] {
            if let Ok(refseq) = RefSeq::from_target(target.clone(), k) {
                sample_ref_kmers.extend(refseq.kmers);
            }
        }

        let adaptive_config = AdaptiveConfig::default();
        let thresholds = estimate_thresholds(db.as_ref(), &sample_ref_kmers, &adaptive_config);

        tracing::info!(
            "adaptive thresholds: tier={}, median_coverage={:.0}, error_rate={:.6}, count={}, ratio={:.4}",
            thresholds.tier,
            thresholds.median_coverage,
            thresholds.error_rate,
            thresholds.count_threshold,
            thresholds.ratio_threshold,
        );

        // CLI --count / --ratio override adaptive values when explicitly provided.
        // clap defaults: count=2, ratio=0.05. We detect "explicit" by checking if
        // the user-supplied value differs from the clap default; if so, they intended
        // an override.
        let final_count = if args.count != 2 {
            tracing::info!("--count {} overrides adaptive count {}", args.count, thresholds.count_threshold);
            args.count
        } else {
            thresholds.count_threshold
        };

        let final_ratio = if (args.ratio - 0.05).abs() > f64::EPSILON {
            tracing::info!("--ratio {} overrides adaptive ratio {}", args.ratio, thresholds.ratio_threshold);
            args.ratio
        } else {
            thresholds.ratio_threshold
        };

        (final_count, final_ratio)
    } else {
        (args.count, args.ratio)
    };

    let walker_config = WalkerConfig {
        count,
        ratio,
        max_stack: args.max_stack,
        max_break: args.max_break,
        max_node: args.max_node,
        adaptive: args.adaptive,
        bidirectional: args.bidirectional,
    };

    if args.bidirectional {
        tracing::info!("bidirectional walking enabled");
    }

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

    // 5b. If --debug-graph is set, ensure the output directory exists.
    if let Some(ref dir) = args.debug_graph {
        std::fs::create_dir_all(dir)
            .with_context(|| format!("creating debug-graph directory: {}", dir.display()))?;
        tracing::info!("writing debug graphs to: {}", dir.display());
    }

    // 5c. Build bootstrap config if enabled.
    let bootstrap_config = if args.bootstrap {
        tracing::info!(
            "bootstrap enabled: {} replicates, {:.0}% CI",
            args.bootstrap_replicates,
            args.bootstrap_confidence * 100.0,
        );
        Some(BootstrapConfig {
            n_replicates: args.bootstrap_replicates,
            confidence_level: args.bootstrap_confidence,
            seed: args.bootstrap_seed,
        })
    } else {
        None
    };

    // 6. Process each target in parallel.
    let debug_dir = args.debug_graph.as_deref();
    let bs_config_ref = bootstrap_config.as_ref();
    let all_calls: Vec<VariantCall> = targets
        .par_iter()
        .flat_map(|target| {
            let result = process_target(
                target,
                db.as_ref(),
                k,
                &walker_config,
                prune_config.as_ref(),
                &db_name,
                debug_dir,
                bs_config_ref,
            );
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

    // 6b. Deduplicate INDEL calls that represent the same biological event.
    let pre_dedup_count = all_calls.len();
    let all_calls = crate::variant::normalize::deduplicate_calls(all_calls);
    let dedup_removed = pre_dedup_count - all_calls.len();
    if dedup_removed > 0 {
        tracing::info!(
            "deduplicated {} equivalent INDEL calls",
            dedup_removed
        );
    }

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
