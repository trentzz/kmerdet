// Output format integration tests.
//
// Tests all output writers (TSV, CSV, VCF, JSON, JSONL) produce valid output
// and that round-trip parsing works correctly.

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
    chrom: &str,
    pos: u64,
    ref_a: &str,
    alt_a: &str,
) -> VariantCall {
    VariantCall {
        sample: sample.to_string(),
        target: target.to_string(),
        variant_type: vtype,
        variant_name: format!("{}:{}/{}:{}", pos, ref_a, alt_a, pos + ref_a.len() as u64),
        rvaf,
        expression: rvaf * 1000.0,
        min_coverage: 42,
        path_score: 42,
        start_kmer_count: 200,
        ref_sequence: "ACGTACGTACGT".to_string(),
        alt_sequence: format!("ACGT{}ACGT", alt_a),
        info: "vs_ref".to_string(),
        chrom: Some(chrom.to_string()),
        pos: Some(pos),
        ref_allele: Some(ref_a.to_string()),
        alt_allele: Some(alt_a.to_string()),
        pvalue: Some(1.5e-10),
        qual: Some(98.24),
        ci_lower: Some(rvaf * 0.8),
        ci_upper: Some(rvaf * 1.2),
    }
}

fn sample_calls() -> Vec<VariantCall> {
    vec![
        make_call("s1", "TP53", VariantType::Substitution, 0.15, "chr17", 7577120, "C", "T"),
        make_call("s1", "NPM1", VariantType::Insertion, 0.25, "chr5", 170837543, "T", "TCTG"),
        make_call("s1", "FLT3", VariantType::Deletion, 0.08, "chr13", 28608250, "GATA", "G"),
        make_call("s1", "KRAS", VariantType::Reference, 1.0, "chr12", 25398284, "G", "G"),
    ]
}

// ---------------------------------------------------------------------------
// TSV output tests
// ---------------------------------------------------------------------------

#[test]
fn test_tsv_correct_columns() {
    let calls = sample_calls();
    let mut buf = Vec::new();
    output::tsv::write(&calls, &mut buf, false).unwrap();
    let output = String::from_utf8(buf).unwrap();

    let lines: Vec<&str> = output.lines().collect();
    assert_eq!(lines.len(), 5, "Header + 4 data lines");

    // Verify header columns
    let header_cols: Vec<&str> = lines[0].split('\t').collect();
    assert_eq!(header_cols.len(), output::DETECT_HEADERS.len());
    for (i, &expected) in output::DETECT_HEADERS.iter().enumerate() {
        assert_eq!(header_cols[i], expected, "Column {} mismatch", i);
    }

    // Verify data is tab-separated with correct column count
    for (i, line) in lines[1..].iter().enumerate() {
        let cols: Vec<&str> = line.split('\t').collect();
        assert_eq!(
            cols.len(),
            output::DETECT_HEADERS.len(),
            "Row {} has wrong number of columns: {} (expected {})",
            i,
            cols.len(),
            output::DETECT_HEADERS.len()
        );
    }
}

#[test]
fn test_tsv_no_header_option() {
    let calls = sample_calls();
    let mut buf = Vec::new();
    output::tsv::write(&calls, &mut buf, true).unwrap();
    let output = String::from_utf8(buf).unwrap();

    let lines: Vec<&str> = output.lines().collect();
    assert_eq!(lines.len(), 4, "No header + 4 data lines");
    // First line should be data, not header
    assert!(lines[0].starts_with("s1"), "First line should be data");
}

#[test]
fn test_tsv_values_correct() {
    let calls = vec![make_call(
        "test_sample", "test_target", VariantType::Substitution, 0.123456,
        "chr1", 42, "A", "T",
    )];
    let mut buf = Vec::new();
    output::tsv::write(&calls, &mut buf, false).unwrap();
    let output = String::from_utf8(buf).unwrap();

    assert!(output.contains("test_sample"));
    assert!(output.contains("test_target"));
    assert!(output.contains("Substitution"));
    assert!(output.contains("0.123456")); // rVAF to 6 decimal places
    assert!(output.contains("chr1"));
    assert!(output.contains("42")); // pos
}

#[test]
fn test_tsv_empty_calls() {
    let calls: Vec<VariantCall> = vec![];
    let mut buf = Vec::new();
    output::tsv::write(&calls, &mut buf, false).unwrap();
    let output = String::from_utf8(buf).unwrap();

    let lines: Vec<&str> = output.lines().collect();
    assert_eq!(lines.len(), 1, "Only header with no data");
    assert!(lines[0].starts_with("sample\t"));
}

// ---------------------------------------------------------------------------
// CSV output tests
// ---------------------------------------------------------------------------

#[test]
fn test_csv_correct_columns() {
    let calls = sample_calls();
    let mut buf = Vec::new();
    output::csv::write(&calls, &mut buf, false).unwrap();
    let output = String::from_utf8(buf).unwrap();

    let lines: Vec<&str> = output.lines().collect();
    assert_eq!(lines.len(), 5, "Header + 4 data lines");

    // CSV header should be comma-separated
    let header_cols: Vec<&str> = lines[0].split(',').collect();
    assert_eq!(header_cols.len(), output::DETECT_HEADERS.len());
    assert_eq!(header_cols[0], "sample");
    assert_eq!(header_cols[1], "target");
}

#[test]
fn test_csv_comma_separated() {
    let calls = sample_calls();
    let mut buf = Vec::new();
    output::csv::write(&calls, &mut buf, false).unwrap();
    let output = String::from_utf8(buf).unwrap();

    // Each data line should have commas
    for line in output.lines().skip(1) {
        assert!(line.contains(','), "CSV data lines should contain commas");
        // Should NOT contain tabs (that would be TSV)
        assert!(!line.contains('\t'), "CSV should not contain tabs");
    }
}

#[test]
fn test_csv_no_header_option() {
    let calls = sample_calls();
    let mut buf = Vec::new();
    output::csv::write(&calls, &mut buf, true).unwrap();
    let output = String::from_utf8(buf).unwrap();

    let lines: Vec<&str> = output.lines().collect();
    assert_eq!(lines.len(), 4, "No header + 4 data lines");
}

// ---------------------------------------------------------------------------
// VCF output tests
// ---------------------------------------------------------------------------

#[test]
fn test_vcf_valid_header() {
    let calls = sample_calls();
    let mut buf = Vec::new();
    output::vcf::write(&calls, &mut buf).unwrap();
    let output = String::from_utf8(buf).unwrap();

    // VCF 4.3 header requirements
    assert!(output.starts_with("##fileformat=VCFv4.3"), "VCF must start with fileformat");
    assert!(output.contains("##source=kmerdet"), "VCF should identify source");
    assert!(output.contains("##INFO=<ID=KVAF"), "VCF should define KVAF info field");
    assert!(output.contains("##INFO=<ID=KCOV"), "VCF should define KCOV info field");
    assert!(output.contains("##INFO=<ID=KEXP"), "VCF should define KEXP info field");
    assert!(output.contains("##INFO=<ID=KTYPE"), "VCF should define KTYPE info field");
    assert!(output.contains("##INFO=<ID=KPV"), "VCF should define KPV info field");
    assert!(
        output.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"),
        "VCF must have column header line"
    );
}

#[test]
fn test_vcf_data_lines() {
    let calls = sample_calls();
    let mut buf = Vec::new();
    output::vcf::write(&calls, &mut buf).unwrap();
    let output = String::from_utf8(buf).unwrap();

    // Data lines start after the #CHROM header line
    let data_lines: Vec<&str> = output
        .lines()
        .filter(|l| !l.starts_with('#'))
        .collect();

    assert_eq!(data_lines.len(), 4, "Should have 4 data lines");

    // Each data line should have 8 tab-separated fields
    for (i, line) in data_lines.iter().enumerate() {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(
            fields.len(),
            8,
            "VCF data line {} should have 8 fields, got {}",
            i,
            fields.len()
        );
        // Field 7 (FILTER) should be PASS
        assert_eq!(fields[6], "PASS", "Filter should be PASS");
        // Field 8 (INFO) should contain KVAF
        assert!(fields[7].contains("KVAF="), "INFO should contain KVAF");
        assert!(fields[7].contains("KCOV="), "INFO should contain KCOV");
        assert!(fields[7].contains("KTYPE="), "INFO should contain KTYPE");
    }

    // Verify specific values for the first call
    let first_fields: Vec<&str> = data_lines[0].split('\t').collect();
    assert_eq!(first_fields[0], "chr17", "CHROM should be chr17");
    assert_eq!(first_fields[1], "7577120", "POS should be 7577120");
    assert_eq!(first_fields[3], "C", "REF should be C");
    assert_eq!(first_fields[4], "T", "ALT should be T");
}

#[test]
fn test_vcf_info_field_format() {
    let calls = vec![make_call(
        "s1", "TP53", VariantType::Substitution, 0.15, "chr17", 100, "A", "T",
    )];
    let mut buf = Vec::new();
    output::vcf::write(&calls, &mut buf).unwrap();
    let output = String::from_utf8(buf).unwrap();

    let data_line = output.lines().find(|l| !l.starts_with('#')).unwrap();
    let info = data_line.split('\t').nth(7).unwrap();

    // INFO field should have semicolon-separated key=value pairs
    assert!(info.contains("KVAF=0.150000"), "KVAF should be formatted correctly");
    assert!(info.contains("KCOV=42"), "KCOV should be present");
    assert!(info.contains("KTYPE=Substitution"), "KTYPE should be present");
    assert!(info.contains("KPV="), "KPV should be present (pvalue was set)");
}

#[test]
fn test_vcf_empty_calls() {
    let calls: Vec<VariantCall> = vec![];
    let mut buf = Vec::new();
    output::vcf::write(&calls, &mut buf).unwrap();
    let output = String::from_utf8(buf).unwrap();

    // Should still have the header
    assert!(output.starts_with("##fileformat=VCFv4.3"));
    assert!(output.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"));
    // No data lines
    let data_lines: Vec<&str> = output
        .lines()
        .filter(|l| !l.starts_with('#'))
        .collect();
    assert!(data_lines.is_empty(), "No data lines for empty calls");
}

// ---------------------------------------------------------------------------
// JSON output tests
// ---------------------------------------------------------------------------

#[test]
fn test_json_valid_array() {
    let calls = sample_calls();
    let mut buf = Vec::new();
    output::json::write_json(&calls, &mut buf).unwrap();
    let output = String::from_utf8(buf).unwrap();

    let parsed: serde_json::Value = serde_json::from_str(&output).unwrap();
    assert!(parsed.is_array(), "JSON output should be an array");
    assert_eq!(parsed.as_array().unwrap().len(), 4);
}

#[test]
fn test_json_contains_all_fields() {
    let calls = vec![make_call(
        "s1", "TP53", VariantType::Substitution, 0.15, "chr17", 100, "A", "T",
    )];
    let mut buf = Vec::new();
    output::json::write_json(&calls, &mut buf).unwrap();
    let output = String::from_utf8(buf).unwrap();

    let parsed: serde_json::Value = serde_json::from_str(&output).unwrap();
    let call = &parsed[0];

    assert_eq!(call["sample"], "s1");
    assert_eq!(call["target"], "TP53");
    assert_eq!(call["variant_type"], "Substitution");
    assert!((call["rvaf"].as_f64().unwrap() - 0.15).abs() < 1e-6);
    assert_eq!(call["min_coverage"], 42);
    assert_eq!(call["chrom"], "chr17");
    assert_eq!(call["pos"], 100);
    assert_eq!(call["ref_allele"], "A");
    assert_eq!(call["alt_allele"], "T");
}

#[test]
fn test_json_empty_calls() {
    let calls: Vec<VariantCall> = vec![];
    let mut buf = Vec::new();
    output::json::write_json(&calls, &mut buf).unwrap();
    let output = String::from_utf8(buf).unwrap();

    let parsed: serde_json::Value = serde_json::from_str(&output).unwrap();
    assert!(parsed.is_array());
    assert_eq!(parsed.as_array().unwrap().len(), 0);
}

#[test]
fn test_json_optional_fields_null_when_absent() {
    let mut call = make_call("s1", "TP53", VariantType::Substitution, 0.15, "chr17", 100, "A", "T");
    call.pvalue = None;
    call.qual = None;
    call.ci_lower = None;
    call.ci_upper = None;
    call.chrom = None;
    call.pos = None;

    let mut buf = Vec::new();
    output::json::write_json(&[call], &mut buf).unwrap();
    let output = String::from_utf8(buf).unwrap();

    let parsed: serde_json::Value = serde_json::from_str(&output).unwrap();
    let obj = &parsed[0];
    assert!(obj["pvalue"].is_null(), "pvalue should be null when absent");
    assert!(obj["qual"].is_null(), "qual should be null when absent");
    assert!(obj["chrom"].is_null(), "chrom should be null when absent");
    assert!(obj["pos"].is_null(), "pos should be null when absent");
}

// ---------------------------------------------------------------------------
// JSONL output tests
// ---------------------------------------------------------------------------

#[test]
fn test_jsonl_one_object_per_line() {
    let calls = sample_calls();
    let mut buf = Vec::new();
    output::json::write_jsonl(&calls, &mut buf).unwrap();
    let output = String::from_utf8(buf).unwrap();

    let lines: Vec<&str> = output.trim().lines().collect();
    assert_eq!(lines.len(), 4, "JSONL should have one line per call");

    // Each line should be valid JSON
    for (i, line) in lines.iter().enumerate() {
        let parsed: Result<serde_json::Value, _> = serde_json::from_str(line);
        assert!(
            parsed.is_ok(),
            "Line {} should be valid JSON: {}",
            i,
            line
        );
        let obj = parsed.unwrap();
        assert!(obj.is_object(), "Each JSONL line should be a JSON object");
    }
}

#[test]
fn test_jsonl_empty_calls() {
    let calls: Vec<VariantCall> = vec![];
    let mut buf = Vec::new();
    output::json::write_jsonl(&calls, &mut buf).unwrap();
    let output = String::from_utf8(buf).unwrap();

    assert!(output.trim().is_empty(), "JSONL with no calls should be empty");
}

#[test]
fn test_jsonl_fields_match_json() {
    // The same call should produce the same field values in JSON and JSONL
    let calls = vec![make_call(
        "s1", "TP53", VariantType::Substitution, 0.15, "chr17", 100, "A", "T",
    )];

    let mut json_buf = Vec::new();
    output::json::write_json(&calls, &mut json_buf).unwrap();
    let json_output = String::from_utf8(json_buf).unwrap();
    let json_parsed: serde_json::Value = serde_json::from_str(&json_output).unwrap();

    let mut jsonl_buf = Vec::new();
    output::json::write_jsonl(&calls, &mut jsonl_buf).unwrap();
    let jsonl_output = String::from_utf8(jsonl_buf).unwrap();
    let jsonl_parsed: serde_json::Value =
        serde_json::from_str(jsonl_output.trim()).unwrap();

    // Compare key fields
    assert_eq!(json_parsed[0]["sample"], jsonl_parsed["sample"]);
    assert_eq!(json_parsed[0]["target"], jsonl_parsed["target"]);
    assert_eq!(json_parsed[0]["rvaf"], jsonl_parsed["rvaf"]);
    assert_eq!(json_parsed[0]["variant_type"], jsonl_parsed["variant_type"]);
}

// ---------------------------------------------------------------------------
// Cross-format consistency tests
// ---------------------------------------------------------------------------

#[test]
fn test_all_formats_same_call_count() {
    let calls = sample_calls();

    // TSV
    let mut tsv_buf = Vec::new();
    output::tsv::write(&calls, &mut tsv_buf, false).unwrap();
    let tsv_data_lines = String::from_utf8(tsv_buf)
        .unwrap()
        .lines()
        .skip(1) // skip header
        .count();

    // CSV
    let mut csv_buf = Vec::new();
    output::csv::write(&calls, &mut csv_buf, false).unwrap();
    let csv_data_lines = String::from_utf8(csv_buf)
        .unwrap()
        .lines()
        .skip(1)
        .count();

    // VCF
    let mut vcf_buf = Vec::new();
    output::vcf::write(&calls, &mut vcf_buf).unwrap();
    let vcf_data_lines = String::from_utf8(vcf_buf)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .count();

    // JSON
    let mut json_buf = Vec::new();
    output::json::write_json(&calls, &mut json_buf).unwrap();
    let json_count = serde_json::from_str::<serde_json::Value>(
        &String::from_utf8(json_buf).unwrap(),
    )
    .unwrap()
    .as_array()
    .unwrap()
    .len();

    // JSONL
    let mut jsonl_buf = Vec::new();
    output::json::write_jsonl(&calls, &mut jsonl_buf).unwrap();
    let jsonl_count = String::from_utf8(jsonl_buf)
        .unwrap()
        .trim()
        .lines()
        .count();

    // All formats should produce the same number of records
    assert_eq!(tsv_data_lines, 4, "TSV data line count");
    assert_eq!(csv_data_lines, 4, "CSV data line count");
    assert_eq!(vcf_data_lines, 4, "VCF data line count");
    assert_eq!(json_count, 4, "JSON array element count");
    assert_eq!(jsonl_count, 4, "JSONL line count");
}

// ---------------------------------------------------------------------------
// TSV round-trip parsing test
// ---------------------------------------------------------------------------

#[test]
fn test_tsv_write_parse_roundtrip() {
    let original = sample_calls();

    // Write to temp file
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("roundtrip.tsv");
    let mut file = std::fs::File::create(&path).unwrap();
    output::tsv::write(&original, &mut file, false).unwrap();
    drop(file);

    // Parse back
    let parsed = output::parse_detection_tsv(&path).unwrap();

    assert_eq!(parsed.len(), original.len());
    for (i, (orig, pars)) in original.iter().zip(parsed.iter()).enumerate() {
        assert_eq!(orig.sample, pars.sample, "sample mismatch at row {}", i);
        assert_eq!(orig.target, pars.target, "target mismatch at row {}", i);
        assert_eq!(
            orig.variant_type, pars.variant_type,
            "variant_type mismatch at row {}",
            i
        );
        assert!(
            (orig.rvaf - pars.rvaf).abs() < 1e-4,
            "rvaf mismatch at row {}: {} vs {}",
            i,
            orig.rvaf,
            pars.rvaf
        );
        assert_eq!(
            orig.min_coverage, pars.min_coverage,
            "min_coverage mismatch at row {}",
            i
        );
        assert_eq!(orig.chrom, pars.chrom, "chrom mismatch at row {}", i);
        assert_eq!(orig.pos, pars.pos, "pos mismatch at row {}", i);
        assert_eq!(
            orig.ref_allele, pars.ref_allele,
            "ref_allele mismatch at row {}",
            i
        );
        assert_eq!(
            orig.alt_allele, pars.alt_allele,
            "alt_allele mismatch at row {}",
            i
        );
    }
}
