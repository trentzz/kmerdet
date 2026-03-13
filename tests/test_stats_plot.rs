// Stats and plot integration tests.
//
// Tests stats computation on detection TSV files and SVG plot generation.

use std::fs;
use std::path::PathBuf;

use kmerdet::output;
use kmerdet::variant::{VariantCall, VariantType};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn make_call(
    sample: &str,
    target: &str,
    vtype: VariantType,
    rvaf: f64,
    expression: f64,
    min_coverage: u64,
) -> VariantCall {
    VariantCall {
        sample: sample.to_string(),
        target: target.to_string(),
        variant_type: vtype,
        variant_name: "100:A/T:100".to_string(),
        rvaf,
        expression,
        min_coverage,
        path_score: min_coverage,
        start_kmer_count: 200,
        ref_sequence: "ACGTACGT".to_string(),
        alt_sequence: "TCGTACGT".to_string(),
        info: "vs_ref".to_string(),
        chrom: Some("chr1".to_string()),
        pos: Some(100),
        ref_allele: Some("A".to_string()),
        alt_allele: Some("T".to_string()),
        pvalue: None,
        qual: None,
        ci_lower: None,
        ci_upper: None,
    }
}

/// Write calls to a TSV file and return the path.
fn write_tsv(dir: &std::path::Path, name: &str, calls: &[VariantCall]) -> PathBuf {
    let path = dir.join(name);
    let mut file = std::fs::File::create(&path).unwrap();
    output::tsv::write(calls, &mut file, false).unwrap();
    path
}

// ---------------------------------------------------------------------------
// Stats tests
// ---------------------------------------------------------------------------

#[test]
fn test_stats_basic_computation() {
    // Create TSV with known values and verify stats output
    let dir = tempfile::tempdir().unwrap();

    let calls = vec![
        make_call("s1", "TP53", VariantType::Substitution, 0.10, 10.0, 50),
        make_call("s1", "BRCA1", VariantType::Reference, 1.0, 100.0, 500),
        make_call("s1", "EGFR", VariantType::Insertion, 0.20, 20.0, 80),
        make_call("s1", "FLT3", VariantType::Itd, 0.30, 30.0, 120),
    ];

    let tsv_path = write_tsv(dir.path(), "detection.tsv", &calls);

    // Run stats via CLI
    let out_path = dir.path().join("stats.json");
    let global = kmerdet::cli::GlobalOptions {
        config: None,
        threads: 0,
        verbose: 0,
        quiet: false,
        output: Some(out_path.clone()),
        format: kmerdet::cli::OutputFormat::Json,
        no_header: false,
    };
    let args = kmerdet::cli::stats::StatsArgs {
        input: vec![tsv_path],
        group_by: None,
    };

    kmerdet::cli::stats::run(args, &global).unwrap();

    // Parse the JSON output
    let content = fs::read_to_string(&out_path).unwrap();
    let stats: serde_json::Value = serde_json::from_str(&content).unwrap();

    // Verify computed statistics
    assert_eq!(stats["total_targets"], 4);
    assert_eq!(stats["detected_variants"], 3); // 3 non-reference
    assert!((stats["detection_rate"].as_f64().unwrap() - 0.75).abs() < 1e-4);

    // VAF mean of non-ref: (0.10 + 0.20 + 0.30) / 3 = 0.2
    let vaf_mean = stats["vaf_mean"].as_f64().unwrap();
    assert!(
        (vaf_mean - 0.2).abs() < 1e-4,
        "VAF mean should be ~0.2, got {}",
        vaf_mean
    );

    // VAF min/max
    let vaf_min = stats["vaf_min"].as_f64().unwrap();
    let vaf_max = stats["vaf_max"].as_f64().unwrap();
    assert!((vaf_min - 0.10).abs() < 1e-4, "VAF min should be 0.10");
    assert!((vaf_max - 0.30).abs() < 1e-4, "VAF max should be 0.30");

    // Variant type counts
    let variant_types = &stats["variant_types"];
    assert_eq!(variant_types["Substitution"], 1);
    assert_eq!(variant_types["Insertion"], 1);
    assert_eq!(variant_types["ITD"], 1);
    // Reference should NOT appear in variant_types
    assert!(
        variant_types.get("Reference").is_none(),
        "Reference should not appear in variant type counts"
    );
}

#[test]
fn test_stats_all_reference_calls() {
    let dir = tempfile::tempdir().unwrap();

    let calls = vec![
        make_call("s1", "TP53", VariantType::Reference, 1.0, 100.0, 500),
        make_call("s1", "BRCA1", VariantType::Reference, 1.0, 80.0, 400),
    ];

    let tsv_path = write_tsv(dir.path(), "all_ref.tsv", &calls);

    let out_path = dir.path().join("stats.json");
    let global = kmerdet::cli::GlobalOptions {
        config: None,
        threads: 0,
        verbose: 0,
        quiet: false,
        output: Some(out_path.clone()),
        format: kmerdet::cli::OutputFormat::Json,
        no_header: false,
    };
    let args = kmerdet::cli::stats::StatsArgs {
        input: vec![tsv_path],
        group_by: None,
    };

    kmerdet::cli::stats::run(args, &global).unwrap();

    let content = fs::read_to_string(&out_path).unwrap();
    let stats: serde_json::Value = serde_json::from_str(&content).unwrap();

    assert_eq!(stats["total_targets"], 2);
    assert_eq!(stats["detected_variants"], 0);
    assert_eq!(stats["detection_rate"], 0.0);
    assert_eq!(stats["vaf_mean"], 0.0);
}

#[test]
fn test_stats_tsv_format() {
    let dir = tempfile::tempdir().unwrap();

    let calls = vec![
        make_call("s1", "TP53", VariantType::Substitution, 0.15, 15.0, 50),
    ];

    let tsv_path = write_tsv(dir.path(), "detection.tsv", &calls);

    let out_path = dir.path().join("stats.tsv");
    let global = kmerdet::cli::GlobalOptions {
        config: None,
        threads: 0,
        verbose: 0,
        quiet: false,
        output: Some(out_path.clone()),
        format: kmerdet::cli::OutputFormat::Tsv,
        no_header: false,
    };
    let args = kmerdet::cli::stats::StatsArgs {
        input: vec![tsv_path],
        group_by: None,
    };

    kmerdet::cli::stats::run(args, &global).unwrap();

    let content = fs::read_to_string(&out_path).unwrap();
    // TSV stats should have metric\tvalue header
    assert!(content.contains("metric\tvalue"), "TSV stats should have header");
    assert!(content.contains("total_targets\t1"));
    assert!(content.contains("detected_variants\t1"));
    assert!(content.contains("vaf_mean\t"));
}

#[test]
fn test_stats_multiple_input_files() {
    let dir = tempfile::tempdir().unwrap();

    let calls1 = vec![
        make_call("s1", "TP53", VariantType::Substitution, 0.10, 10.0, 50),
    ];
    let calls2 = vec![
        make_call("s2", "EGFR", VariantType::Insertion, 0.20, 20.0, 80),
    ];

    let f1 = write_tsv(dir.path(), "s1.tsv", &calls1);
    let f2 = write_tsv(dir.path(), "s2.tsv", &calls2);

    let out_path = dir.path().join("combined_stats.json");
    let global = kmerdet::cli::GlobalOptions {
        config: None,
        threads: 0,
        verbose: 0,
        quiet: false,
        output: Some(out_path.clone()),
        format: kmerdet::cli::OutputFormat::Json,
        no_header: false,
    };
    let args = kmerdet::cli::stats::StatsArgs {
        input: vec![f1, f2],
        group_by: None,
    };

    kmerdet::cli::stats::run(args, &global).unwrap();

    let content = fs::read_to_string(&out_path).unwrap();
    let stats: serde_json::Value = serde_json::from_str(&content).unwrap();

    // Combined stats should reflect both files
    assert_eq!(stats["total_targets"], 2);
    assert_eq!(stats["detected_variants"], 2);
}

#[test]
fn test_stats_group_by_sample() {
    let dir = tempfile::tempdir().unwrap();

    let calls = vec![
        make_call("s1", "TP53", VariantType::Substitution, 0.10, 10.0, 50),
        make_call("s1", "BRCA1", VariantType::Reference, 1.0, 100.0, 500),
        make_call("s2", "TP53", VariantType::Substitution, 0.20, 20.0, 80),
        make_call("s2", "EGFR", VariantType::Insertion, 0.30, 30.0, 120),
    ];

    let tsv_path = write_tsv(dir.path(), "multi_sample.tsv", &calls);

    let out_path = dir.path().join("grouped_stats.tsv");
    let global = kmerdet::cli::GlobalOptions {
        config: None,
        threads: 0,
        verbose: 0,
        quiet: false,
        output: Some(out_path.clone()),
        format: kmerdet::cli::OutputFormat::Tsv,
        no_header: false,
    };
    let args = kmerdet::cli::stats::StatsArgs {
        input: vec![tsv_path],
        group_by: Some("sample".to_string()),
    };

    kmerdet::cli::stats::run(args, &global).unwrap();

    let content = fs::read_to_string(&out_path).unwrap();
    // Should have header + 2 group rows (s1, s2)
    let lines: Vec<&str> = content.lines().collect();
    assert!(lines.len() >= 3, "Should have header + 2 sample groups");
    assert!(content.contains("s1"), "Should include s1 group");
    assert!(content.contains("s2"), "Should include s2 group");
}

#[test]
fn test_stats_empty_input_fails() {
    let global = kmerdet::cli::GlobalOptions {
        config: None,
        threads: 0,
        verbose: 0,
        quiet: false,
        output: None,
        format: kmerdet::cli::OutputFormat::Json,
        no_header: false,
    };
    let args = kmerdet::cli::stats::StatsArgs {
        input: vec![],
        group_by: None,
    };

    let result = kmerdet::cli::stats::run(args, &global);
    assert!(result.is_err(), "Empty input should fail");
}

// ---------------------------------------------------------------------------
// Plot SVG generation tests
// ---------------------------------------------------------------------------

#[test]
fn test_plot_vaf_histogram_generates_valid_svg() {
    let dir = tempfile::tempdir().unwrap();

    let calls = vec![
        make_call("s1", "TP53", VariantType::Substitution, 0.10, 10.0, 50),
        make_call("s1", "EGFR", VariantType::Insertion, 0.20, 20.0, 80),
        make_call("s1", "FLT3", VariantType::Itd, 0.30, 30.0, 120),
        make_call("s1", "BRCA1", VariantType::Reference, 1.0, 100.0, 500),
    ];

    let tsv_path = write_tsv(dir.path(), "detection.tsv", &calls);
    let svg_path = dir.path().join("vaf_histogram.svg");

    let args = kmerdet::cli::plot::PlotArgs {
        input: tsv_path,
        chart: kmerdet::cli::plot::ChartType::VafHistogram,
        output: Some(svg_path.clone()),
    };

    let global = kmerdet::cli::GlobalOptions {
        config: None,
        threads: 0,
        verbose: 0,
        quiet: false,
        output: None,
        format: kmerdet::cli::OutputFormat::Tsv,
        no_header: false,
    };

    kmerdet::cli::plot::run(args, &global).unwrap();

    let svg = fs::read_to_string(&svg_path).unwrap();
    assert!(svg.starts_with("<svg"), "Output should be valid SVG");
    assert!(svg.contains("</svg>"), "SVG should be properly closed");
    assert!(svg.contains("VAF Distribution"), "Should contain title");
    assert!(svg.contains("rect"), "Should contain bar rectangles");
}

#[test]
fn test_plot_type_distribution_generates_valid_svg() {
    let dir = tempfile::tempdir().unwrap();

    let calls = vec![
        make_call("s1", "TP53", VariantType::Substitution, 0.10, 10.0, 50),
        make_call("s1", "BRCA1", VariantType::Substitution, 0.15, 15.0, 60),
        make_call("s1", "EGFR", VariantType::Insertion, 0.20, 20.0, 80),
        make_call("s1", "FLT3", VariantType::Itd, 0.30, 30.0, 120),
        make_call("s1", "NPM1", VariantType::Deletion, 0.08, 8.0, 30),
    ];

    let tsv_path = write_tsv(dir.path(), "detection.tsv", &calls);
    let svg_path = dir.path().join("type_dist.svg");

    let args = kmerdet::cli::plot::PlotArgs {
        input: tsv_path,
        chart: kmerdet::cli::plot::ChartType::TypeDistribution,
        output: Some(svg_path.clone()),
    };

    let global = kmerdet::cli::GlobalOptions {
        config: None,
        threads: 0,
        verbose: 0,
        quiet: false,
        output: None,
        format: kmerdet::cli::OutputFormat::Tsv,
        no_header: false,
    };

    kmerdet::cli::plot::run(args, &global).unwrap();

    let svg = fs::read_to_string(&svg_path).unwrap();
    assert!(svg.starts_with("<svg"));
    assert!(svg.contains("Variant Type Distribution"));
    assert!(svg.contains("</svg>"));
    // Should mention variant types
    assert!(svg.contains("Substitution"));
    assert!(svg.contains("Insertion"));
}

#[test]
fn test_plot_detection_bar_generates_valid_svg() {
    let dir = tempfile::tempdir().unwrap();

    let calls = vec![
        make_call("s1", "targetA", VariantType::Reference, 1.0, 100.0, 500),
        make_call("s1", "targetA", VariantType::Substitution, 0.15, 15.0, 50),
        make_call("s1", "targetB", VariantType::Insertion, 0.08, 8.0, 30),
        make_call("s1", "targetB", VariantType::Reference, 1.0, 80.0, 400),
    ];

    let tsv_path = write_tsv(dir.path(), "detection.tsv", &calls);
    let svg_path = dir.path().join("detection_bar.svg");

    let args = kmerdet::cli::plot::PlotArgs {
        input: tsv_path,
        chart: kmerdet::cli::plot::ChartType::DetectionBar,
        output: Some(svg_path.clone()),
    };

    let global = kmerdet::cli::GlobalOptions {
        config: None,
        threads: 0,
        verbose: 0,
        quiet: false,
        output: None,
        format: kmerdet::cli::OutputFormat::Tsv,
        no_header: false,
    };

    kmerdet::cli::plot::run(args, &global).unwrap();

    let svg = fs::read_to_string(&svg_path).unwrap();
    assert!(svg.starts_with("<svg"));
    assert!(svg.contains("Detection by Target"));
    assert!(svg.contains("</svg>"));
}

#[test]
fn test_plot_summary_pie_generates_valid_svg() {
    let dir = tempfile::tempdir().unwrap();

    let calls = vec![
        make_call("s1", "t1", VariantType::Substitution, 0.15, 15.0, 50),
        make_call("s1", "t2", VariantType::Reference, 1.0, 100.0, 500),
        make_call("s1", "t3", VariantType::Reference, 1.0, 80.0, 400),
    ];

    let tsv_path = write_tsv(dir.path(), "detection.tsv", &calls);
    let svg_path = dir.path().join("summary_pie.svg");

    let args = kmerdet::cli::plot::PlotArgs {
        input: tsv_path,
        chart: kmerdet::cli::plot::ChartType::SummaryPie,
        output: Some(svg_path.clone()),
    };

    let global = kmerdet::cli::GlobalOptions {
        config: None,
        threads: 0,
        verbose: 0,
        quiet: false,
        output: None,
        format: kmerdet::cli::OutputFormat::Tsv,
        no_header: false,
    };

    kmerdet::cli::plot::run(args, &global).unwrap();

    let svg = fs::read_to_string(&svg_path).unwrap();
    assert!(svg.starts_with("<svg"));
    assert!(svg.contains("Detection Summary"));
    assert!(svg.contains("</svg>"));
    // Pie chart should have path or circle elements
    assert!(
        svg.contains("path") || svg.contains("circle"),
        "Pie chart should contain path or circle elements"
    );
}

#[test]
fn test_plot_with_only_reference_calls_fails() {
    let dir = tempfile::tempdir().unwrap();

    let calls = vec![
        make_call("s1", "t1", VariantType::Reference, 1.0, 100.0, 500),
    ];

    let tsv_path = write_tsv(dir.path(), "all_ref.tsv", &calls);
    let svg_path = dir.path().join("should_fail.svg");

    let args = kmerdet::cli::plot::PlotArgs {
        input: tsv_path,
        chart: kmerdet::cli::plot::ChartType::VafHistogram,
        output: Some(svg_path),
    };

    let global = kmerdet::cli::GlobalOptions {
        config: None,
        threads: 0,
        verbose: 0,
        quiet: false,
        output: None,
        format: kmerdet::cli::OutputFormat::Tsv,
        no_header: false,
    };

    // VAF histogram with only reference calls should fail
    let result = kmerdet::cli::plot::run(args, &global);
    assert!(
        result.is_err(),
        "VAF histogram with only reference calls should fail"
    );
}
