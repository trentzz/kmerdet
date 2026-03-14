use std::path::PathBuf;

use anyhow::{Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;

use crate::graph::pruning::PruneConfig;
use crate::jellyfish::KmerDatabase;
use crate::sequence::path::KmerPath;
use crate::sequence::target::{RefSeq, Target};
use crate::trace::{self, DetectionTrace};
use crate::variant::bootstrap::BootstrapConfig;
use crate::variant::classifier::Classification;
use crate::variant::cluster;
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

    /// Write per-target detection traces to DIR as JSON files
    #[arg(long, value_name = "DIR")]
    pub trace: Option<PathBuf>,

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

    /// Run detection at multiple k-mer lengths and merge with consensus voting.
    ///
    /// Specify comma-separated k values (e.g., `--multi-k 21,31,41`).
    /// Shorter k values improve INDEL sensitivity; longer k values improve SNV
    /// specificity. Results are merged using weighted consensus voting.
    #[arg(long, value_delimiter = ',', num_args = 1..)]
    pub multi_k: Option<Vec<u8>>,

    /// Write diagnostic report to this directory.
    #[arg(long, value_name = "DIR")]
    pub report_dir: Option<PathBuf>,

    /// Diagnostic report verbosity level [default: summary]
    #[arg(long, value_enum, default_value = "summary")]
    pub report_level: crate::report::ReportLevel,

    /// Sensitivity preset: ultra, high, standard, strict.
    ///
    /// Sets coordinated defaults for count, ratio, min-qual, min-rvaf, and
    /// min-coverage. Explicit CLI flags override the preset values.
    /// Can also be configured in the [sensitivity] section of the config file.
    #[arg(long, value_name = "PRESET")]
    pub sensitivity: Option<String>,

    /// Minimum QUAL (Phred-scaled p-value) for reporting a variant.
    /// Calls with QUAL below this are excluded from output.
    #[arg(long)]
    pub min_qual: Option<f64>,

    /// Minimum rVAF for reporting a variant.
    /// Calls with rVAF below this are excluded from output.
    #[arg(long)]
    pub min_rvaf: Option<f64>,

    /// Minimum k-mer coverage for reporting a variant.
    /// Calls with min_coverage below this are excluded from output.
    #[arg(long)]
    pub min_cov: Option<u64>,
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
///
/// When `trace_enabled` is true, also builds and returns a `DetectionTrace` capturing
/// the intermediate results for debugging/audit purposes.
fn process_target(
    target: &Target,
    db: &dyn KmerDatabase,
    k: u8,
    walker_config: &WalkerConfig,
    prune_config: Option<&PruneConfig>,
    sample: &str,
    debug_graph_dir: Option<&std::path::Path>,
    bootstrap_config: Option<&BootstrapConfig>,
    trace_enabled: bool,
) -> Result<(Vec<VariantCall>, Option<DetectionTrace>)> {
    // a. Decompose target into reference k-mers.
    let refseq = RefSeq::from_target(target.clone(), k)?;

    // b. Walk the k-mer graph starting from reference k-mers.
    let walk_result = if walker_config.bidirectional {
        crate::walker::walk_bidirectional(db, &refseq.kmers, walker_config)
    } else {
        crate::walker::walk(db, &refseq.kmers, walker_config)
    };

    // Build walking trace if tracing is enabled.
    let walking_trace = if trace_enabled {
        Some(trace::build_walking_trace(&walk_result, &refseq.kmers, db, walker_config))
    } else {
        None
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

    // Build graph trace if tracing is enabled.
    let graph_trace = if trace_enabled {
        Some(trace::build_graph_trace(&graph))
    } else {
        None
    };

    // d. Find alternative paths ranked by coverage (highest-coverage alt first).
    let ranked_paths = crate::graph::pathfind::find_alternative_paths_ranked(&graph, db);

    // d'. Optionally write a debug DOT file for this target.
    if let Some(dir) = debug_graph_dir {
        let safe_name: String = target
            .name
            .chars()
            .map(|c| if c.is_alphanumeric() || c == '-' || c == '_' { c } else { '_' })
            .collect();
        let dot_path = dir.join(format!("{}.dot", safe_name));

        let paths_for_dot: Vec<KmerPath> = ranked_paths.iter().map(|(p, _)| p.clone()).collect();
        let path_refs: Vec<&crate::sequence::path::KmerPath> = paths_for_dot.iter().collect();
        let highlight = if paths_for_dot.is_empty() {
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

    if ranked_paths.is_empty() {
        // No paths at all (disconnected graph) -- produce a Reference call anyway.
        let calls = vec![make_reference_call(sample, target, db, &refseq)];

        let dt = if trace_enabled {
            Some(DetectionTrace {
                target: target.name.clone(),
                walking: walking_trace.unwrap(),
                graph: graph_trace.unwrap(),
                pathfinding: trace::build_pathfinding_trace(&[]),
                classifications: vec![],
                quantification: None,
                outcome: "no_paths".to_string(),
            })
        } else {
            None
        };

        return Ok((calls, dt));
    }

    // Extract paths (preserving coverage-ranked order) for quantification.
    let paths: Vec<KmerPath> = ranked_paths.iter().map(|(p, _)| p.clone()).collect();

    // The first path is always the reference path.
    let ref_path = &paths[0];

    if paths.len() == 1 {
        // Only reference path: no variants detected.
        let calls = vec![make_reference_call(sample, target, db, &refseq)];

        let dt = if trace_enabled {
            Some(DetectionTrace {
                target: target.name.clone(),
                walking: walking_trace.unwrap(),
                graph: graph_trace.unwrap(),
                pathfinding: trace::build_pathfinding_trace(&paths),
                classifications: vec![],
                quantification: None,
                outcome: "reference_only".to_string(),
            })
        } else {
            None
        };

        return Ok((calls, dt));
    }

    // Log coverage scores for debugging.
    for (i, (_, score)) in ranked_paths.iter().enumerate() {
        if i == 0 {
            tracing::debug!(
                target: "detect",
                "ref path coverage score: {}",
                score
            );
        } else {
            tracing::debug!(
                target: "detect",
                "alt path {} coverage score: {}",
                i,
                score
            );
        }
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
    //    Alternative paths are already sorted by coverage (highest first),
    //    so the primary variant call (index 1) is the most confident.
    let mut calls: Vec<VariantCall> = Vec::new();
    let mut classification_traces: Vec<trace::ClassificationTrace> = Vec::new();

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
        ranked_paths[0].1,
        db,
    ));

    // Process each alternative path (indices 1..).
    for (i, alt_path) in paths.iter().enumerate().skip(1) {
        let classification =
            crate::variant::classifier::classify(ref_path, alt_path, k);

        // Build classification trace if tracing is enabled.
        if trace_enabled {
            classification_traces.push(trace::build_classification_trace(&classification, i));
        }

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
            ranked_paths[i].1,
            db,
        );
        call.pvalue = Some(pvalue);
        call.qual = Some(qual);

        calls.push(call);
    }

    // Build the full detection trace if enabled.
    let dt = if trace_enabled {
        let quant_trace = trace::build_quantification_trace(&quant, &paths);
        Some(DetectionTrace {
            target: target.name.clone(),
            walking: walking_trace.unwrap(),
            graph: graph_trace.unwrap(),
            pathfinding: trace::build_pathfinding_trace(&paths),
            classifications: classification_traces,
            quantification: Some(quant_trace),
            outcome: "variant_detected".to_string(),
        })
    } else {
        None
    };

    Ok((calls, dt))
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
    path_score: u64,
    db: &dyn KmerDatabase,
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
        start_kmer_count: alt_path.kmers.first().map_or(0, |k| db.query(k)),
        ref_sequence: ref_path.to_sequence(),
        alt_sequence: alt_path.to_sequence(),
        info: "vs_ref".to_string(),
        path_score,
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
        path_score: min_count,
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

impl DetectArgs {
    /// Apply values from a config file for any args that were not explicitly
    /// set on the command line.
    pub(crate) fn apply_config(
        &mut self,
        cfg: &crate::config::Config,
        matches: &clap::ArgMatches,
    ) {
        let is_set = |id| super::is_subcommand_arg_set(matches, id);

        // db: only apply if not provided on CLI (it's required, so this is
        // only useful when the config supplies it and the user omits it --
        // but clap will error first for required args. We handle it for
        // completeness in case we later make db optional.)
        if self.targets.is_empty() {
            if !cfg.detect.targets.is_empty() {
                self.targets = cfg.detect.targets.clone();
            }
        }

        if !is_set("count") {
            if let Some(v) = cfg.detect.count {
                self.count = v;
            }
        }
        if !is_set("ratio") {
            if let Some(v) = cfg.detect.ratio {
                self.ratio = v;
            }
        }
        if !is_set("max_stack") {
            if let Some(v) = cfg.detect.max_stack {
                self.max_stack = v;
            }
        }
        if !is_set("max_break") {
            if let Some(v) = cfg.detect.max_break {
                self.max_break = v;
            }
        }
        if !is_set("max_node") {
            if let Some(v) = cfg.detect.max_node {
                self.max_node = v;
            }
        }
        if !is_set("cluster") {
            if let Some(v) = cfg.detect.cluster {
                self.cluster = v;
            }
        }

        // Apply [sensitivity] section: preset first, then individual overrides.
        if !is_set("sensitivity") {
            if cfg.sensitivity.preset.is_some() {
                self.sensitivity = cfg.sensitivity.preset.clone();
            }
        }
        if !is_set("min_qual") {
            if let Some(v) = cfg.sensitivity.min_qual {
                self.min_qual = Some(v);
            }
        }
        if !is_set("min_rvaf") {
            if let Some(v) = cfg.sensitivity.min_rvaf {
                self.min_rvaf = Some(v);
            }
        }
        if !is_set("min_cov") {
            if let Some(v) = cfg.sensitivity.min_coverage {
                self.min_cov = Some(v);
            }
        }
        // Sensitivity section can also set count/ratio if not already set by
        // [detect] or CLI.
        if !is_set("count") && cfg.detect.count.is_none() {
            if let Some(v) = cfg.sensitivity.count {
                self.count = v;
            }
        }
        if !is_set("ratio") && cfg.detect.ratio.is_none() {
            if let Some(v) = cfg.sensitivity.ratio {
                self.ratio = v;
            }
        }
        if !is_set("adaptive") {
            if let Some(true) = cfg.sensitivity.adaptive {
                self.adaptive = true;
            }
        }
        if !is_set("bidirectional") {
            if let Some(true) = cfg.sensitivity.bidirectional {
                self.bidirectional = true;
            }
        }
        if !is_set("bootstrap") {
            if let Some(true) = cfg.sensitivity.bootstrap {
                self.bootstrap = true;
            }
        }
    }

    /// Apply a sensitivity preset, overriding count/ratio and setting
    /// min_qual/min_rvaf/min_cov defaults if not already explicitly set.
    pub(crate) fn apply_sensitivity_preset(&mut self, matches: &clap::ArgMatches) -> Result<()> {
        let preset_name = match &self.sensitivity {
            Some(name) => name.clone(),
            None => return Ok(()),
        };

        let preset = crate::config::SensitivityPresetValues::from_name(&preset_name)
            .ok_or_else(|| {
                anyhow::anyhow!(
                    "unknown sensitivity preset '{}'; valid presets: {}",
                    preset_name,
                    crate::config::SensitivityPresetValues::preset_names().join(", ")
                )
            })?;

        let is_set = |id| super::is_subcommand_arg_set(matches, id);

        // Preset overrides count/ratio unless explicitly set on CLI.
        if !is_set("count") {
            self.count = preset.count;
        }
        if !is_set("ratio") {
            self.ratio = preset.ratio;
        }
        if self.min_qual.is_none() && !is_set("min_qual") {
            self.min_qual = Some(preset.min_qual);
        }
        if self.min_rvaf.is_none() && !is_set("min_rvaf") {
            self.min_rvaf = Some(preset.min_rvaf);
        }
        if self.min_cov.is_none() && !is_set("min_cov") {
            self.min_cov = Some(preset.min_coverage);
        }

        tracing::info!(
            "sensitivity preset '{}': count={}, ratio={}, min_qual={}, min_rvaf={}, min_cov={}",
            preset_name,
            self.count,
            self.ratio,
            self.min_qual.unwrap_or(0.0),
            self.min_rvaf.unwrap_or(0.0),
            self.min_cov.unwrap_or(0),
        );

        Ok(())
    }
}

pub fn run(
    mut args: DetectArgs,
    global: &super::GlobalOptions,
    cfg: Option<&crate::config::Config>,
    matches: &clap::ArgMatches,
) -> Result<()> {
    if let Some(cfg) = cfg {
        args.apply_config(cfg, matches);
    }

    // Apply sensitivity preset (after config, so CLI > config > preset > defaults).
    args.apply_sensitivity_preset(matches)?;

    // Check for multi-k mode first.
    if let Some(k_values) = args.multi_k.clone() {
        return run_multi_k(args, global, k_values);
    }

    run_single_k(args, global)
}

/// Run detection at a single k-mer length (the default mode).
fn run_single_k(args: DetectArgs, global: &super::GlobalOptions) -> Result<()> {
    let pipeline_start = std::time::Instant::now();

    // 1. Open jellyfish database.
    let db = open_db(&args.db)?;
    let db_name = sample_name(&args.db);

    // 2. Determine k-mer length (CLI override > DB header).
    let k = resolve_kmer_length(args.kmer_length, &args.db)?;
    tracing::info!("using k-mer length: {}", k);

    // 3. Load all target FASTA files.
    let targets = load_all_targets(&args.targets)?;
    let n_targets = targets.len();

    // 4. Build walker configuration from args.
    let walker_config = build_walker_config(&args, db.as_ref(), &targets, k)?;

    if args.bidirectional {
        tracing::info!("bidirectional walking enabled");
    }

    // 4b. Build prune configuration (unless --no-prune).
    let prune_config = build_prune_config(&args);

    // 5. Run detection pipeline.
    let (all_calls, traces) = run_detection_pipeline(
        &args,
        db.as_ref(),
        &db_name,
        k,
        &targets,
        &walker_config,
        prune_config.as_ref(),
    )?;

    // 6. Apply sensitivity thresholds (min_qual, min_rvaf, min_cov).
    let all_calls = apply_sensitivity_thresholds(all_calls, &args);

    // 7. Write diagnostic report if --report-dir is set.
    if let Some(ref report_dir) = args.report_dir {
        let total_time_ms = pipeline_start.elapsed().as_millis() as u64;
        if let Err(e) = write_report(
            report_dir,
            args.report_level,
            &all_calls,
            &traces,
            n_targets,
            total_time_ms,
            k,
            &args,
        ) {
            tracing::warn!("failed to write diagnostic report: {}", e);
        } else {
            tracing::info!("wrote diagnostic report to {}", report_dir.display());
        }
    }

    // 8. Write output.
    write_output(&all_calls, global)?;

    Ok(())
}

/// Run detection at multiple k-mer lengths and merge with consensus voting.
fn run_multi_k(
    args: DetectArgs,
    global: &super::GlobalOptions,
    k_values: Vec<u8>,
) -> Result<()> {
    anyhow::ensure!(
        !k_values.is_empty(),
        "--multi-k requires at least one k value"
    );

    // Validate k values.
    for &k in &k_values {
        anyhow::ensure!(k >= 5, "k-mer length {} is too small (minimum 5)", k);
        anyhow::ensure!(k <= 127, "k-mer length {} exceeds maximum (127)", k);
        anyhow::ensure!(k % 2 == 1, "k-mer length {} must be odd for canonical form", k);
    }

    let total_k_values = k_values.len();
    tracing::info!(
        "multi-k detection enabled: k={:?} ({} values)",
        k_values,
        total_k_values
    );

    // 1. Open jellyfish database.
    let db = open_db(&args.db)?;
    let db_name = sample_name(&args.db);

    // 2. Load all target FASTA files (shared across all k values).
    let targets = load_all_targets(&args.targets)?;

    // 3. Build prune configuration (shared across all k values).
    let prune_config = build_prune_config(&args);

    // 4. Run detection for each k value independently.
    let mut calls_per_k: Vec<(u8, Vec<VariantCall>)> = Vec::new();

    for &k in &k_values {
        tracing::info!("--- running detection at k={} ---", k);

        // Build walker config for this k value.
        let walker_config = build_walker_config(&args, db.as_ref(), &targets, k)?;

        // Run detection at this k value.
        let (calls, _traces) = run_detection_pipeline(
            &args,
            db.as_ref(),
            &db_name,
            k,
            &targets,
            &walker_config,
            prune_config.as_ref(),
        )?;

        tracing::info!("k={}: detected {} calls", k, calls.len());
        calls_per_k.push((k, calls));
    }

    // 5. Merge results with consensus voting.
    let consensus_calls =
        crate::variant::consensus::merge_multi_k(&calls_per_k, total_k_values);

    tracing::info!(
        "consensus merge: {} unique variants from {} total calls across {} k values",
        consensus_calls.len(),
        calls_per_k.iter().map(|(_, c)| c.len()).sum::<usize>(),
        total_k_values,
    );

    // Log tier distribution.
    let tier_counts = consensus_calls.iter().fold([0usize; 4], |mut acc, c| {
        match c.tier {
            crate::variant::consensus::ConfidenceTier::Tier1 => acc[0] += 1,
            crate::variant::consensus::ConfidenceTier::Tier2 => acc[1] += 1,
            crate::variant::consensus::ConfidenceTier::Tier3 => acc[2] += 1,
            crate::variant::consensus::ConfidenceTier::Tier4 => acc[3] += 1,
        }
        acc
    });
    tracing::info!(
        "tier distribution: Tier1={}, Tier2={}, Tier3={}, Tier4={}",
        tier_counts[0],
        tier_counts[1],
        tier_counts[2],
        tier_counts[3],
    );

    // 6. Convert consensus calls to VariantCalls for output.
    let all_calls: Vec<VariantCall> = consensus_calls
        .iter()
        .map(crate::variant::consensus::consensus_to_variant_call)
        .collect();

    // 7. Write output.
    write_output(&all_calls, global)?;

    Ok(())
}

/// Load all target FASTA files from the specified paths.
fn load_all_targets(paths: &[PathBuf]) -> Result<Vec<Target>> {
    let mut targets: Vec<Target> = Vec::new();
    for path in paths {
        let loaded = crate::sequence::target::load_targets(path)
            .with_context(|| format!("loading targets from {}", path.display()))?;
        targets.extend(loaded);
    }
    tracing::info!("loaded {} targets", targets.len());

    if targets.is_empty() {
        anyhow::bail!("no targets loaded; check --targets paths");
    }

    Ok(targets)
}

/// Build the WalkerConfig from CLI args, potentially using adaptive thresholds.
fn build_walker_config(
    args: &DetectArgs,
    db: &dyn KmerDatabase,
    targets: &[Target],
    k: u8,
) -> Result<WalkerConfig> {
    let (count, ratio) = if args.adaptive {
        let mut sample_ref_kmers: Vec<String> = Vec::new();
        let max_targets_to_sample = 10.min(targets.len());
        for target in &targets[..max_targets_to_sample] {
            if let Ok(refseq) = RefSeq::from_target(target.clone(), k) {
                sample_ref_kmers.extend(refseq.kmers);
            }
        }

        let adaptive_config = AdaptiveConfig::default();
        let thresholds = estimate_thresholds(db, &sample_ref_kmers, &adaptive_config);

        tracing::info!(
            "adaptive thresholds (k={}): tier={}, median_coverage={:.0}, error_rate={:.6}, count={}, ratio={:.4}",
            k,
            thresholds.tier,
            thresholds.median_coverage,
            thresholds.error_rate,
            thresholds.count_threshold,
            thresholds.ratio_threshold,
        );

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

    Ok(WalkerConfig {
        count,
        ratio,
        max_stack: args.max_stack,
        max_break: args.max_break,
        max_node: args.max_node,
        adaptive: args.adaptive,
        bidirectional: args.bidirectional,
    })
}

/// Build the PruneConfig from CLI args.
fn build_prune_config(args: &DetectArgs) -> Option<PruneConfig> {
    if args.no_prune {
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
    }
}

/// Run the detection pipeline for a single k value: walk, graph, classify, quantify.
/// Returns variant calls and optional detection traces (when tracing is enabled).
fn run_detection_pipeline(
    args: &DetectArgs,
    db: &dyn KmerDatabase,
    db_name: &str,
    k: u8,
    targets: &[Target],
    walker_config: &WalkerConfig,
    prune_config: Option<&PruneConfig>,
) -> Result<(Vec<VariantCall>, Vec<DetectionTrace>)> {
    // Set up progress bar.
    let pb = ProgressBar::new(targets.len() as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template(
                &format!(
                    "{{spinner:.green}} [k={}] [{{elapsed_precise}}] [{{bar:40.cyan/blue}}] {{pos}}/{{len}} targets ({{eta}})",
                    k
                ),
            )
            .unwrap()
            .progress_chars("=>-"),
    );

    // If --debug-graph is set, ensure the output directory exists.
    if let Some(ref dir) = args.debug_graph {
        std::fs::create_dir_all(dir)
            .with_context(|| format!("creating debug-graph directory: {}", dir.display()))?;
    }

    // Enable tracing if --trace or --report-dir is set.
    let trace_enabled = args.trace.is_some() || args.report_dir.is_some();
    if let Some(ref dir) = args.trace {
        std::fs::create_dir_all(dir)
            .with_context(|| format!("creating trace directory: {}", dir.display()))?;
    }

    // Build bootstrap config if enabled.
    let bootstrap_config = if args.bootstrap {
        Some(BootstrapConfig {
            n_replicates: args.bootstrap_replicates,
            confidence_level: args.bootstrap_confidence,
            seed: args.bootstrap_seed,
        })
    } else {
        None
    };

    // Process each target in parallel, collecting both calls and optional traces.
    let debug_dir = args.debug_graph.as_deref();
    let bs_config_ref = bootstrap_config.as_ref();
    let results: Vec<(Vec<VariantCall>, Option<DetectionTrace>)> = targets
        .par_iter()
        .filter_map(|target| {
            let result = process_target(
                target,
                db,
                k,
                walker_config,
                prune_config,
                db_name,
                debug_dir,
                bs_config_ref,
                trace_enabled,
            );
            pb.inc(1);
            match result {
                Ok(pair) => Some(pair),
                Err(e) => {
                    tracing::warn!("target '{}' failed at k={}: {}", target.name, k, e);
                    None
                }
            }
        })
        .collect();

    pb.finish_with_message("done");

    // Separate calls and traces.
    let mut all_calls: Vec<VariantCall> = Vec::new();
    let mut all_traces: Vec<DetectionTrace> = Vec::new();

    for (calls, dt) in results {
        all_calls.extend(calls);
        if let Some(t) = dt {
            all_traces.push(t);
        }
    }

    // Write per-target traces and summary if --trace is enabled.
    if let Some(ref trace_dir) = args.trace {
        for dt in &all_traces {
            if let Err(e) = trace::write_trace_target(dt, trace_dir) {
                tracing::warn!("failed to write trace for target '{}': {}", dt.target, e);
            }
        }
        if let Err(e) = trace::write_trace_summary(&all_traces, trace_dir) {
            tracing::warn!("failed to write trace summary: {}", e);
        } else {
            tracing::info!(
                "wrote detection traces for {} targets to {}",
                all_traces.len(),
                trace_dir.display(),
            );
        }
    }

    // Apply overlapping mutation clustering if --cluster is enabled.
    if args.cluster {
        let clusters = cluster::cluster_variants(&all_calls);
        let multi_clusters: Vec<&cluster::VariantCluster> =
            clusters.iter().filter(|c| c.calls.len() > 1).collect();

        if !multi_clusters.is_empty() {
            tracing::info!(
                "found {} overlapping mutation cluster(s):",
                multi_clusters.len(),
            );
            for (i, c) in multi_clusters.iter().enumerate() {
                tracing::info!(
                    "  cluster {}: {} calls spanning positions [{}, {}]",
                    i + 1,
                    c.calls.len(),
                    c.start,
                    c.end,
                );
                for call in &c.calls {
                    tracing::info!(
                        "    - {} {} (rVAF={:.4})",
                        call.variant_type,
                        call.variant_name,
                        call.rvaf,
                    );
                }
            }

            // Annotate clustered calls: update the info field to record cluster membership.
            // Build a lookup: variant_name -> cluster index for calls that belong to a
            // multi-variant cluster.
            let mut cluster_membership: std::collections::HashMap<String, usize> =
                std::collections::HashMap::new();
            for (i, c) in multi_clusters.iter().enumerate() {
                for call in &c.calls {
                    cluster_membership.insert(call.variant_name.clone(), i + 1);
                }
            }

            for call in &mut all_calls {
                if let Some(cluster_id) = cluster_membership.get(&call.variant_name) {
                    if !call.info.is_empty() {
                        call.info.push_str(&format!(";cluster={}", cluster_id));
                    } else {
                        call.info = format!("cluster={}", cluster_id);
                    }
                }
            }
        }
    }

    // Deduplicate INDEL calls.
    let pre_dedup_count = all_calls.len();
    let all_calls = crate::variant::normalize::deduplicate_calls(all_calls);
    let dedup_removed = pre_dedup_count - all_calls.len();
    if dedup_removed > 0 {
        tracing::info!(
            "k={}: deduplicated {} equivalent INDEL calls",
            k,
            dedup_removed,
        );
    }

    tracing::info!(
        "k={}: detected {} variant calls across {} targets",
        k,
        all_calls.len(),
        targets.len(),
    );

    Ok((all_calls, all_traces))
}

/// Apply sensitivity thresholds (min_qual, min_rvaf, min_cov) to filter
/// variant calls. Reference calls are always kept.
fn apply_sensitivity_thresholds(calls: Vec<VariantCall>, args: &DetectArgs) -> Vec<VariantCall> {
    let min_qual = args.min_qual.unwrap_or(0.0);
    let min_rvaf = args.min_rvaf.unwrap_or(0.0);
    let min_cov = args.min_cov.unwrap_or(0);

    // Skip filtering if all thresholds are at their most permissive.
    if min_qual <= 0.0 && min_rvaf <= 0.0 && min_cov == 0 {
        return calls;
    }

    let pre_count = calls
        .iter()
        .filter(|c| c.variant_type != VariantType::Reference)
        .count();

    let filtered: Vec<VariantCall> = calls
        .into_iter()
        .filter(|call| {
            // Always keep reference calls.
            if call.variant_type == VariantType::Reference {
                return true;
            }
            // Apply QUAL threshold.
            if min_qual > 0.0 {
                if let Some(qual) = call.qual {
                    if qual < min_qual {
                        return false;
                    }
                } else {
                    // No QUAL computed; exclude if threshold is set.
                    return false;
                }
            }
            // Apply rVAF threshold.
            if min_rvaf > 0.0 && call.rvaf < min_rvaf {
                return false;
            }
            // Apply coverage threshold.
            if min_cov > 0 && call.min_coverage < min_cov {
                return false;
            }
            true
        })
        .collect();

    let post_count = filtered
        .iter()
        .filter(|c| c.variant_type != VariantType::Reference)
        .count();
    let removed = pre_count - post_count;
    if removed > 0 {
        tracing::info!(
            "sensitivity thresholds removed {} variant calls (min_qual={}, min_rvaf={}, min_cov={})",
            removed,
            min_qual,
            min_rvaf,
            min_cov,
        );
    }

    filtered
}

/// Write variant calls to the output destination.
fn write_output(all_calls: &[VariantCall], global: &super::GlobalOptions) -> Result<()> {
    match &global.output {
        Some(path) => {
            if global.format == super::OutputFormat::Xlsx {
                crate::output::write_calls_xlsx(all_calls, path)?;
                return Ok(());
            }
            let mut file = std::fs::File::create(path)
                .with_context(|| format!("creating output file: {}", path.display()))?;
            crate::output::write_calls(all_calls, &global.format, &mut file, global.no_header)?;
        }
        None => {
            let mut stdout = std::io::stdout().lock();
            crate::output::write_calls(
                all_calls,
                &global.format,
                &mut stdout,
                global.no_header,
            )?;
        }
    }
    Ok(())
}

/// Write diagnostic report from detection results, using trace data when available.
fn write_report(
    report_dir: &std::path::Path,
    report_level: crate::report::ReportLevel,
    all_calls: &[VariantCall],
    traces: &[DetectionTrace],
    n_targets: usize,
    total_time_ms: u64,
    k: u8,
    args: &DetectArgs,
) -> Result<()> {
    use crate::report::*;

    let writer = ReportWriter::new(report_dir.to_path_buf(), report_level)?;

    // Build a lookup from target name to trace for enriching reports.
    let trace_by_target: std::collections::HashMap<&str, &DetectionTrace> = traces
        .iter()
        .map(|t| (t.target.as_str(), t))
        .collect();

    // Group calls by target to build per-target summaries.
    let mut calls_by_target: indexmap::IndexMap<String, Vec<&VariantCall>> =
        indexmap::IndexMap::new();
    for call in all_calls {
        calls_by_target
            .entry(call.target.clone())
            .or_default()
            .push(call);
    }

    // Count how many targets failed (processed but produced no calls at all).
    let targets_with_calls = calls_by_target.len();
    let targets_failed = n_targets.saturating_sub(targets_with_calls);

    // Build detection summary entries and count stats.
    let mut summary_entries: Vec<DetectionSummaryEntry> = Vec::new();
    let mut total_variants = 0usize;
    let mut targets_with_variants = 0usize;
    let mut targets_reference_only = 0usize;

    for (target_name, calls) in &calls_by_target {
        let variant_calls: Vec<&&VariantCall> = calls
            .iter()
            .filter(|c| c.variant_type != VariantType::Reference)
            .collect();
        let n_variants = variant_calls.len();
        total_variants += n_variants;

        let result_type = if let Some(vc) = variant_calls.first() {
            vc.variant_type.to_string()
        } else {
            "Reference".to_string()
        };

        let rvaf = variant_calls
            .iter()
            .map(|c| c.rvaf)
            .reduce(f64::max)
            .unwrap_or(1.0);

        let coverage = calls.first().map(|c| c.min_coverage).unwrap_or(0);
        let n_paths = calls.len();

        if n_variants > 0 {
            targets_with_variants += 1;
        } else {
            targets_reference_only += 1;
        }

        summary_entries.push(DetectionSummaryEntry {
            target: target_name.clone(),
            result_type: result_type.clone(),
            rvaf,
            coverage,
            n_paths,
            n_variants,
        });

        // Write per-target result at standard+ level.
        if report_level >= ReportLevel::Standard {
            let variants: Vec<TargetVariant> = variant_calls
                .iter()
                .map(|c| TargetVariant {
                    variant_type: c.variant_type.to_string(),
                    variant_name: c.variant_name.clone(),
                    rvaf: c.rvaf,
                    min_coverage: c.min_coverage,
                    ref_allele: c.ref_allele.clone().unwrap_or_default(),
                    alt_allele: c.alt_allele.clone().unwrap_or_default(),
                })
                .collect();

            let decision_text = build_decision_text(&result_type, &variant_calls);

            // Use trace data to populate walking/graph/path stats.
            let trace = trace_by_target.get(target_name.as_str());
            let target_result = TargetResult {
                target: target_name.clone(),
                result_type,
                variants,
                n_ref_kmers: trace.map_or(0, |t| t.walking.reference_kmers),
                n_walked_kmers: trace.map_or(0, |t| t.walking.total_nodes),
                n_graph_nodes: trace.map_or(0, |t| t.graph.total_nodes),
                n_paths,
                decision_text,
            };

            writer.write_target_result(target_name, &target_result)?;
        }
    }

    // Write run summary (always written).
    let run_summary = RunSummary {
        targets_processed: n_targets,
        variants_found: total_variants,
        targets_with_variants,
        targets_reference_only,
        targets_failed,
        total_time_ms,
        parameters: RunParameters {
            kmer_length: k,
            count: args.count,
            ratio: args.ratio,
            adaptive: args.adaptive,
            bidirectional: args.bidirectional,
            prune_enabled: !args.no_prune,
        },
    };
    writer.write_run_summary(&run_summary)?;

    // Write detection summary TSV (always written).
    writer.write_detection_summary(&summary_entries)?;

    Ok(())
}

/// Build a human-readable decision text for a target.
fn build_decision_text(result_type: &str, variant_calls: &[&&VariantCall]) -> String {
    if variant_calls.is_empty() {
        return "Reference only -- no variant branches found".to_string();
    }

    let mut lines: Vec<String> = Vec::new();
    for vc in variant_calls {
        let pos_info = if let Some(pos) = vc.pos {
            format!(" at pos {}", pos)
        } else if !vc.variant_name.is_empty() {
            format!(" {}", vc.variant_name)
        } else {
            String::new()
        };
        lines.push(format!(
            "Detected {}{}, rVAF={:.4}, coverage={}",
            vc.variant_type, pos_info, vc.rvaf, vc.min_coverage
        ));
    }

    if lines.len() == 1 {
        lines.into_iter().next().unwrap()
    } else {
        format!(
            "Multiple variants detected ({} {}):\n{}",
            variant_calls.len(),
            result_type,
            lines.join("\n")
        )
    }
}
