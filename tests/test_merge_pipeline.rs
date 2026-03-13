// Merge pipeline integration tests.
//
// Tests the merge subcommand by creating multiple detection TSV files,
// merging them, and verifying the output.

use std::fs::{self, File};
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
) -> VariantCall {
    VariantCall {
        sample: sample.to_string(),
        target: target.to_string(),
        variant_type: vtype,
        variant_name: format!("100:A/T:100"),
        rvaf,
        expression: rvaf * 100.0,
        min_coverage: 50,
        path_score: 50,
        start_kmer_count: 100,
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

/// Write a set of VariantCalls to a TSV file and return the path.
fn write_tsv(dir: &std::path::Path, name: &str, calls: &[VariantCall]) -> PathBuf {
    let path = dir.join(name);
    let mut file = File::create(&path).unwrap();
    output::tsv::write(calls, &mut file, false).unwrap();
    path
}

/// Read a TSV file line-by-line and return all lines.
fn read_lines(path: &std::path::Path) -> Vec<String> {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .map(|s| s.to_string())
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[test]
fn test_merge_two_files_via_tsv() {
    let dir = tempfile::tempdir().unwrap();

    let calls1 = vec![
        make_call("s1", "TP53", VariantType::Substitution, 0.15),
        make_call("s1", "BRCA1", VariantType::Reference, 1.0),
    ];
    let calls2 = vec![
        make_call("s2", "TP53", VariantType::Substitution, 0.20),
        make_call("s2", "EGFR", VariantType::Insertion, 0.08),
    ];

    let f1 = write_tsv(dir.path(), "sample1.tsv", &calls1);
    let f2 = write_tsv(dir.path(), "sample2.tsv", &calls2);

    // Merge via CLI merge subcommand logic
    let out_path = dir.path().join("merged.tsv");
    let global = kmerdet::cli::GlobalOptions {
        config: None,
        threads: 0,
        verbose: 0,
        quiet: false,
        output: Some(out_path.clone()),
        format: kmerdet::cli::OutputFormat::Tsv,
        no_header: false,
    };
    let args = kmerdet::cli::merge::MergeArgs {
        input: vec![f1, f2],
        deduplicate: false,
    };

    kmerdet::cli::merge::run(args, &global).unwrap();

    let lines = read_lines(&out_path);
    // Header + 4 data lines
    assert_eq!(lines.len(), 5, "Merged file should have header + 4 rows");
    assert!(lines[0].starts_with("sample\t"), "First line should be header");
    assert!(lines[1].contains("s1"), "Row 1 should be from sample 1");
    assert!(lines[2].contains("s1"), "Row 2 should be from sample 1");
    assert!(lines[3].contains("s2"), "Row 3 should be from sample 2");
    assert!(lines[4].contains("s2"), "Row 4 should be from sample 2");
}

#[test]
fn test_merge_three_files() {
    let dir = tempfile::tempdir().unwrap();

    let f1 = write_tsv(
        dir.path(),
        "a.tsv",
        &[make_call("s1", "TP53", VariantType::Substitution, 0.10)],
    );
    let f2 = write_tsv(
        dir.path(),
        "b.tsv",
        &[make_call("s2", "BRCA1", VariantType::Deletion, 0.20)],
    );
    let f3 = write_tsv(
        dir.path(),
        "c.tsv",
        &[make_call("s3", "FLT3", VariantType::Itd, 0.30)],
    );

    let out_path = dir.path().join("merged.tsv");
    let global = kmerdet::cli::GlobalOptions {
        config: None,
        threads: 0,
        verbose: 0,
        quiet: false,
        output: Some(out_path.clone()),
        format: kmerdet::cli::OutputFormat::Tsv,
        no_header: false,
    };
    let args = kmerdet::cli::merge::MergeArgs {
        input: vec![f1, f2, f3],
        deduplicate: false,
    };

    kmerdet::cli::merge::run(args, &global).unwrap();

    let lines = read_lines(&out_path);
    assert_eq!(lines.len(), 4, "Header + 3 data rows");
}

#[test]
fn test_merge_deduplication() {
    let dir = tempfile::tempdir().unwrap();

    // Same sample/target/variant_name in both files (different rVAFs)
    let f1 = write_tsv(
        dir.path(),
        "a.tsv",
        &[make_call("s1", "TP53", VariantType::Substitution, 0.10)],
    );
    let f2 = write_tsv(
        dir.path(),
        "b.tsv",
        &[make_call("s1", "TP53", VariantType::Substitution, 0.15)],
    );

    let out_path = dir.path().join("deduped.tsv");
    let global = kmerdet::cli::GlobalOptions {
        config: None,
        threads: 0,
        verbose: 0,
        quiet: false,
        output: Some(out_path.clone()),
        format: kmerdet::cli::OutputFormat::Tsv,
        no_header: false,
    };
    let args = kmerdet::cli::merge::MergeArgs {
        input: vec![f1, f2],
        deduplicate: true,
    };

    kmerdet::cli::merge::run(args, &global).unwrap();

    let lines = read_lines(&out_path);
    assert_eq!(
        lines.len(),
        2,
        "Deduplication should keep header + 1 row (duplicate removed)"
    );
    // First occurrence wins
    assert!(
        lines[1].contains("0.100000"),
        "First occurrence (rVAF=0.10) should be kept"
    );
}

#[test]
fn test_merge_no_dedup_keeps_all() {
    let dir = tempfile::tempdir().unwrap();

    let f1 = write_tsv(
        dir.path(),
        "a.tsv",
        &[make_call("s1", "TP53", VariantType::Substitution, 0.10)],
    );
    let f2 = write_tsv(
        dir.path(),
        "b.tsv",
        &[make_call("s1", "TP53", VariantType::Substitution, 0.15)],
    );

    let out_path = dir.path().join("not_deduped.tsv");
    let global = kmerdet::cli::GlobalOptions {
        config: None,
        threads: 0,
        verbose: 0,
        quiet: false,
        output: Some(out_path.clone()),
        format: kmerdet::cli::OutputFormat::Tsv,
        no_header: false,
    };
    let args = kmerdet::cli::merge::MergeArgs {
        input: vec![f1, f2],
        deduplicate: false,
    };

    kmerdet::cli::merge::run(args, &global).unwrap();

    let lines = read_lines(&out_path);
    assert_eq!(
        lines.len(),
        3,
        "Without dedup, both rows should be kept (header + 2)"
    );
}

#[test]
fn test_merge_header_mismatch_fails() {
    let dir = tempfile::tempdir().unwrap();

    // Write files with different headers
    let f1 = dir.path().join("a.tsv");
    let f2 = dir.path().join("b.tsv");
    fs::write(&f1, "col_a\tcol_b\nv1\tv2\n").unwrap();
    fs::write(&f2, "col_x\tcol_y\nv3\tv4\n").unwrap();

    let global = kmerdet::cli::GlobalOptions {
        config: None,
        threads: 0,
        verbose: 0,
        quiet: false,
        output: None,
        format: kmerdet::cli::OutputFormat::Tsv,
        no_header: false,
    };
    let args = kmerdet::cli::merge::MergeArgs {
        input: vec![f1, f2],
        deduplicate: false,
    };

    let result = kmerdet::cli::merge::run(args, &global);
    assert!(result.is_err(), "Mismatched headers should cause an error");
    assert!(
        result.unwrap_err().to_string().contains("header mismatch"),
        "Error should mention header mismatch"
    );
}

#[test]
fn test_merge_empty_input_fails() {
    let global = kmerdet::cli::GlobalOptions {
        config: None,
        threads: 0,
        verbose: 0,
        quiet: false,
        output: None,
        format: kmerdet::cli::OutputFormat::Tsv,
        no_header: false,
    };
    let args = kmerdet::cli::merge::MergeArgs {
        input: vec![],
        deduplicate: false,
    };

    let result = kmerdet::cli::merge::run(args, &global);
    assert!(result.is_err(), "Empty input should fail");
}

#[test]
fn test_merge_preserves_all_columns() {
    let dir = tempfile::tempdir().unwrap();

    let calls = vec![VariantCall {
        sample: "test_sample".to_string(),
        target: "NPM1_4ins".to_string(),
        variant_type: VariantType::Insertion,
        variant_name: "41:T/TCTG:41".to_string(),
        rvaf: 0.25,
        expression: 150.0,
        min_coverage: 42,
        path_score: 42,
        start_kmer_count: 200,
        ref_sequence: "ACGTACGT".to_string(),
        alt_sequence: "ACGTCTGACGT".to_string(),
        info: "vs_ref".to_string(),
        chrom: Some("chr5".to_string()),
        pos: Some(170837543),
        ref_allele: Some("T".to_string()),
        alt_allele: Some("TCTG".to_string()),
        pvalue: Some(1.5e-12),
        qual: Some(118.24),
        ci_lower: Some(0.20),
        ci_upper: Some(0.30),
    }];

    let f1 = write_tsv(dir.path(), "full.tsv", &calls);

    let out_path = dir.path().join("merged.tsv");
    let global = kmerdet::cli::GlobalOptions {
        config: None,
        threads: 0,
        verbose: 0,
        quiet: false,
        output: Some(out_path.clone()),
        format: kmerdet::cli::OutputFormat::Tsv,
        no_header: false,
    };
    let args = kmerdet::cli::merge::MergeArgs {
        input: vec![f1],
        deduplicate: false,
    };

    kmerdet::cli::merge::run(args, &global).unwrap();

    // Parse the merged output and verify all fields
    let parsed = output::parse_detection_tsv(&out_path).unwrap();
    assert_eq!(parsed.len(), 1);
    assert_eq!(parsed[0].sample, "test_sample");
    assert_eq!(parsed[0].target, "NPM1_4ins");
    assert_eq!(parsed[0].variant_type, VariantType::Insertion);
    assert!((parsed[0].rvaf - 0.25).abs() < 1e-4);
    assert_eq!(parsed[0].min_coverage, 42);
    assert_eq!(parsed[0].chrom.as_deref(), Some("chr5"));
    assert_eq!(parsed[0].pos, Some(170837543));
}
