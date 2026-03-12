use kmerdet::output;
use kmerdet::variant::{VariantCall, VariantType};

fn sample_call() -> VariantCall {
    VariantCall {
        sample: "test_sample".to_string(),
        target: "NPM1_4ins".to_string(),
        variant_type: VariantType::Insertion,
        variant_name: "41:T/TCTG:41".to_string(),
        rvaf: 0.25,
        expression: 150.0,
        min_coverage: 42,
        start_kmer_count: 200,
        ref_sequence: "ACGTACGT".to_string(),
        alt_sequence: "ACGTCTGACGT".to_string(),
        info: "vs_ref".to_string(),
        chrom: Some("chr5".to_string()),
        pos: Some(170837543),
        ref_allele: Some("T".to_string()),
        alt_allele: Some("TCTG".to_string()),
        ci_lower: None,
        ci_upper: None,
    }
}

#[test]
fn test_tsv_output_with_header() {
    let calls = vec![sample_call()];
    let mut buf = Vec::new();
    output::tsv::write(&calls, &mut buf, false).unwrap();
    let output = String::from_utf8(buf).unwrap();

    assert!(output.starts_with("sample\t"));
    assert!(output.contains("test_sample"));
    assert!(output.contains("Insertion"));
    assert!(output.contains("0.250000"));
}

#[test]
fn test_tsv_output_no_header() {
    let calls = vec![sample_call()];
    let mut buf = Vec::new();
    output::tsv::write(&calls, &mut buf, true).unwrap();
    let output = String::from_utf8(buf).unwrap();

    assert!(!output.starts_with("sample\t"));
    assert!(output.starts_with("test_sample"));
}

#[test]
fn test_json_output() {
    let calls = vec![sample_call()];
    let mut buf = Vec::new();
    output::json::write_json(&calls, &mut buf).unwrap();
    let output = String::from_utf8(buf).unwrap();

    let parsed: serde_json::Value = serde_json::from_str(&output).unwrap();
    assert!(parsed.is_array());
    assert_eq!(parsed.as_array().unwrap().len(), 1);
}

#[test]
fn test_jsonl_output() {
    let calls = vec![sample_call(), sample_call()];
    let mut buf = Vec::new();
    output::json::write_jsonl(&calls, &mut buf).unwrap();
    let output = String::from_utf8(buf).unwrap();

    let lines: Vec<&str> = output.trim().lines().collect();
    assert_eq!(lines.len(), 2);
    // Each line should be valid JSON
    for line in lines {
        serde_json::from_str::<serde_json::Value>(line).unwrap();
    }
}

#[test]
fn test_vcf_output_header() {
    let calls = vec![sample_call()];
    let mut buf = Vec::new();
    output::vcf::write(&calls, &mut buf).unwrap();
    let output = String::from_utf8(buf).unwrap();

    assert!(output.starts_with("##fileformat=VCFv4.3"));
    assert!(output.contains("##INFO=<ID=KVAF"));
    assert!(output.contains("#CHROM\tPOS\tID\tREF\tALT"));
    assert!(output.contains("chr5\t170837543"));
}

#[test]
fn test_csv_output() {
    let calls = vec![sample_call()];
    let mut buf = Vec::new();
    output::csv::write(&calls, &mut buf, false).unwrap();
    let output = String::from_utf8(buf).unwrap();

    assert!(output.contains("sample,"));
    assert!(output.contains("test_sample"));
}
