//! CLI subcommand for ground truth benchmarking.
//!
//! Compares kmerdet detection results against a known ground truth file
//! and reports sensitivity, specificity, precision, F1-score, and other
//! performance metrics.

use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

use crate::benchmark::{self, GroundTruthVariant};
use crate::variant::VariantCall;

#[derive(clap::Args, Debug)]
pub struct BenchmarkArgs {
    /// Detection results file (TSV from detect)
    #[arg(short, long)]
    pub results: PathBuf,

    /// Ground truth file (TSV: chrom, pos, ref, alt, type, true_vaf, [category])
    #[arg(short, long)]
    pub truth: PathBuf,

    /// VAF bin boundaries (comma-separated)
    #[arg(long, default_value = "0,0.001,0.01,0.1,1.0")]
    pub vaf_bins: String,

    /// VAF threshold sweep values (comma-separated)
    #[arg(long)]
    pub sweep_vaf: Option<String>,
}

pub fn run(args: BenchmarkArgs, global: &super::GlobalOptions) -> Result<()> {
    // 1. Parse detection results from TSV
    let calls = parse_detection_results(&args.results)?;

    // 2. Parse ground truth from TSV
    let truth = parse_ground_truth(&args.truth)?;

    // 3. Parse VAF bin boundaries
    let vaf_bins = parse_f64_list(&args.vaf_bins)
        .context("parsing --vaf-bins")?;
    anyhow::ensure!(
        vaf_bins.len() >= 2,
        "--vaf-bins must have at least 2 boundaries"
    );

    // 4. Parse optional sweep thresholds
    let sweep_thresholds = match &args.sweep_vaf {
        Some(s) => Some(parse_f64_list(s).context("parsing --sweep-vaf")?),
        None => None,
    };

    // 5. Run benchmark
    let report = benchmark::run_benchmark(
        &calls,
        &truth,
        &vaf_bins,
        sweep_thresholds.as_deref(),
    );

    // 6. Write output
    let mut output: Box<dyn Write> = match &global.output {
        Some(path) => Box::new(
            File::create(path)
                .with_context(|| format!("creating output file: {}", path.display()))?,
        ),
        None => Box::new(std::io::stdout().lock()),
    };

    match global.format {
        super::OutputFormat::Json => {
            benchmark::write_benchmark_json(&report, &mut output)?;
        }
        _ => {
            benchmark::write_benchmark_text(&report, &mut output)?;
        }
    }

    Ok(())
}

/// Parse detection results from a TSV file into VariantCall structs.
///
/// Uses the same column-index-based parsing as `filter.rs` and `stats.rs`.
fn parse_detection_results(path: &Path) -> Result<Vec<VariantCall>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("opening detection results: {}", path.display()))?;

    let mut calls = Vec::new();
    for (line_no, record) in reader.records().enumerate() {
        let record =
            record.with_context(|| format!("reading record at line {}", line_no + 2))?;
        let call = VariantCall {
            sample: record.get(0).unwrap_or("").to_string(),
            target: record.get(1).unwrap_or("").to_string(),
            variant_type: record
                .get(2)
                .unwrap_or("Reference")
                .parse()
                .with_context(|| {
                    format!(
                        "parsing variant_type at line {}: {:?}",
                        line_no + 2,
                        record.get(2)
                    )
                })?,
            variant_name: record.get(3).unwrap_or("").to_string(),
            rvaf: record.get(4).unwrap_or("0").parse().unwrap_or(0.0),
            expression: record.get(5).unwrap_or("0").parse().unwrap_or(0.0),
            min_coverage: record.get(6).unwrap_or("0").parse().unwrap_or(0),
            start_kmer_count: record.get(7).unwrap_or("0").parse().unwrap_or(0),
            ref_sequence: record.get(8).unwrap_or("").to_string(),
            alt_sequence: record.get(9).unwrap_or("").to_string(),
            info: record.get(10).unwrap_or("").to_string(),
            chrom: record
                .get(11)
                .and_then(|s| if s.is_empty() { None } else { Some(s.to_string()) }),
            pos: record.get(12).and_then(|s| s.parse().ok()),
            ref_allele: record
                .get(13)
                .and_then(|s| if s.is_empty() { None } else { Some(s.to_string()) }),
            alt_allele: record
                .get(14)
                .and_then(|s| if s.is_empty() { None } else { Some(s.to_string()) }),
            pvalue: record.get(15).and_then(|s| s.parse().ok()),
            qual: record.get(16).and_then(|s| s.parse().ok()),
            ci_lower: record.get(17).and_then(|s| s.parse().ok()),
            ci_upper: record.get(18).and_then(|s| s.parse().ok()),
        };
        calls.push(call);
    }
    Ok(calls)
}

/// Parse ground truth variants from a TSV file.
///
/// Expected format (with header):
/// `chrom\tpos\tref\talt\ttype\ttrue_vaf[\tcategory]`
///
/// The `category` column is optional.
fn parse_ground_truth(path: &Path) -> Result<Vec<GroundTruthVariant>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("opening ground truth file: {}", path.display()))?;

    let mut variants = Vec::new();
    for (line_no, record) in reader.records().enumerate() {
        let record =
            record.with_context(|| format!("reading truth record at line {}", line_no + 2))?;

        let chrom = record
            .get(0)
            .unwrap_or("")
            .to_string();
        let pos: u64 = record
            .get(1)
            .unwrap_or("0")
            .parse()
            .with_context(|| format!("parsing pos at line {}", line_no + 2))?;
        let ref_allele = record
            .get(2)
            .unwrap_or("")
            .to_string();
        let alt_allele = record
            .get(3)
            .unwrap_or("")
            .to_string();
        let variant_type = record
            .get(4)
            .unwrap_or("")
            .to_string();
        let true_vaf: f64 = record
            .get(5)
            .unwrap_or("0")
            .parse()
            .with_context(|| format!("parsing true_vaf at line {}", line_no + 2))?;
        let category = record
            .get(6)
            .and_then(|s| if s.is_empty() { None } else { Some(s.to_string()) });

        variants.push(GroundTruthVariant {
            chrom,
            pos,
            ref_allele,
            alt_allele,
            variant_type,
            true_vaf,
            category,
        });
    }
    Ok(variants)
}

/// Parse a comma-separated list of f64 values.
fn parse_f64_list(s: &str) -> Result<Vec<f64>> {
    s.split(',')
        .map(|v| {
            v.trim()
                .parse::<f64>()
                .with_context(|| format!("invalid number: '{}'", v.trim()))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Create a temporary TSV file with detection results.
    fn write_detection_tsv(records: &[&str]) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        // Header
        writeln!(
            f,
            "sample\ttarget\tvariant_type\tvariant_name\trvaf\texpression\t\
             min_coverage\tstart_kmer_count\tref_sequence\talt_sequence\t\
             info\tchrom\tpos\tref_allele\talt_allele"
        )
        .unwrap();
        for record in records {
            writeln!(f, "{}", record).unwrap();
        }
        f.flush().unwrap();
        f
    }

    /// Create a temporary TSV file with ground truth.
    fn write_truth_tsv(records: &[&str]) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, "chrom\tpos\tref\talt\ttype\ttrue_vaf\tcategory").unwrap();
        for record in records {
            writeln!(f, "{}", record).unwrap();
        }
        f.flush().unwrap();
        f
    }

    #[test]
    fn test_parse_detection_results() {
        let f = write_detection_tsv(&[
            "s1\ttarget1\tSubstitution\tA>T\t0.05\t1.5\t10\t20\tACGT\tTCGT\tvs_ref\tchr1\t100\tA\tT",
        ]);
        let calls = parse_detection_results(f.path()).unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].sample, "s1");
        assert_eq!(calls[0].chrom, Some("chr1".to_string()));
        assert_eq!(calls[0].pos, Some(100));
        assert!((calls[0].rvaf - 0.05).abs() < 1e-10);
    }

    #[test]
    fn test_parse_detection_results_empty() {
        let f = write_detection_tsv(&[]);
        let calls = parse_detection_results(f.path()).unwrap();
        assert!(calls.is_empty());
    }

    #[test]
    fn test_parse_detection_results_missing_optional_fields() {
        let f = write_detection_tsv(&[
            "s1\ttarget1\tReference\tref\t0.0\t0.0\t5\t10\tACGT\tACGT\tvs_ref\t\t\t\t",
        ]);
        let calls = parse_detection_results(f.path()).unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].chrom, None);
        assert_eq!(calls[0].pos, None);
        assert_eq!(calls[0].ref_allele, None);
        assert_eq!(calls[0].alt_allele, None);
    }

    #[test]
    fn test_parse_ground_truth() {
        let f = write_truth_tsv(&[
            "chr1\t100\tA\tT\tSubstitution\t0.05\tSNV",
            "chr2\t200\tG\tGACC\tInsertion\t0.10\t",
        ]);
        let truth = parse_ground_truth(f.path()).unwrap();
        assert_eq!(truth.len(), 2);
        assert_eq!(truth[0].chrom, "chr1");
        assert_eq!(truth[0].pos, 100);
        assert_eq!(truth[0].ref_allele, "A");
        assert_eq!(truth[0].alt_allele, "T");
        assert_eq!(truth[0].variant_type, "Substitution");
        assert!((truth[0].true_vaf - 0.05).abs() < 1e-10);
        assert_eq!(truth[0].category, Some("SNV".to_string()));
        assert_eq!(truth[1].category, None);
    }

    #[test]
    fn test_parse_ground_truth_without_category_column() {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(f, "chrom\tpos\tref\talt\ttype\ttrue_vaf").unwrap();
        writeln!(f, "chr1\t100\tA\tT\tSubstitution\t0.05").unwrap();
        f.flush().unwrap();

        let truth = parse_ground_truth(f.path()).unwrap();
        assert_eq!(truth.len(), 1);
        assert_eq!(truth[0].category, None);
    }

    #[test]
    fn test_parse_ground_truth_absent_variant() {
        let f = write_truth_tsv(&["chr1\t100\tA\tT\tSubstitution\t0.0\t"]);
        let truth = parse_ground_truth(f.path()).unwrap();
        assert_eq!(truth[0].true_vaf, 0.0);
    }

    #[test]
    fn test_parse_f64_list() {
        let result = parse_f64_list("0,0.001,0.01,0.1,1.0").unwrap();
        assert_eq!(result, vec![0.0, 0.001, 0.01, 0.1, 1.0]);
    }

    #[test]
    fn test_parse_f64_list_with_spaces() {
        let result = parse_f64_list("0.0, 0.5, 1.0").unwrap();
        assert_eq!(result, vec![0.0, 0.5, 1.0]);
    }

    #[test]
    fn test_parse_f64_list_invalid() {
        let result = parse_f64_list("0.0,abc,1.0");
        assert!(result.is_err());
    }

    #[test]
    fn test_parse_f64_list_single() {
        let result = parse_f64_list("0.5").unwrap();
        assert_eq!(result, vec![0.5]);
    }

    #[test]
    fn test_end_to_end_benchmark_text() {
        let det_file = write_detection_tsv(&[
            "s1\ttarget1\tSubstitution\tA>T\t0.05\t1.5\t10\t20\tACGT\tTCGT\tvs_ref\tchr1\t100\tA\tT",
            "s1\ttarget2\tInsertion\tins\t0.10\t2.0\t15\t25\tGGGG\tGGGACCG\tvs_ref\tchr2\t200\tG\tGACC",
        ]);
        let truth_file = write_truth_tsv(&[
            "chr1\t100\tA\tT\tSubstitution\t0.05\tSNV",
            "chr2\t200\tG\tGACC\tInsertion\t0.10\tINS",
            "chr3\t300\tC\tT\tSubstitution\t0.03\tSNV",
        ]);

        let calls = parse_detection_results(det_file.path()).unwrap();
        let truth = parse_ground_truth(truth_file.path()).unwrap();
        let vaf_bins = vec![0.0, 0.01, 0.1, 1.0];

        let report = benchmark::run_benchmark(&calls, &truth, &vaf_bins, None);

        assert_eq!(report.confusion.tp, 2);
        assert_eq!(report.confusion.fn_count, 1);
        assert_eq!(report.confusion.fp, 0);

        let mut buf = Vec::new();
        benchmark::write_benchmark_text(&report, &mut buf).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("Sensitivity"));
        assert!(output.contains("Per-Type Breakdown:"));
    }

    #[test]
    fn test_end_to_end_benchmark_json() {
        let det_file = write_detection_tsv(&[
            "s1\ttarget1\tSubstitution\tA>T\t0.05\t1.5\t10\t20\tACGT\tTCGT\tvs_ref\tchr1\t100\tA\tT",
        ]);
        let truth_file = write_truth_tsv(&[
            "chr1\t100\tA\tT\tSubstitution\t0.05\t",
        ]);

        let calls = parse_detection_results(det_file.path()).unwrap();
        let truth = parse_ground_truth(truth_file.path()).unwrap();

        let report = benchmark::run_benchmark(&calls, &truth, &[0.0, 1.0], None);

        let mut buf = Vec::new();
        benchmark::write_benchmark_json(&report, &mut buf).unwrap();
        let output = String::from_utf8(buf).unwrap();

        let parsed: serde_json::Value = serde_json::from_str(&output).unwrap();
        assert_eq!(parsed["confusion"]["tp"], 1);
        assert_eq!(parsed["confusion"]["fp"], 0);
    }

    #[test]
    fn test_end_to_end_with_sweep() {
        let det_file = write_detection_tsv(&[
            "s1\ttarget1\tSubstitution\tA>T\t0.05\t1.5\t10\t20\tACGT\tTCGT\tvs_ref\tchr1\t100\tA\tT",
            "s1\ttarget2\tSubstitution\tG>C\t0.001\t0.5\t3\t5\tGGGG\tCGGG\tvs_ref\tchr2\t200\tG\tC",
        ]);
        let truth_file = write_truth_tsv(&[
            "chr1\t100\tA\tT\tSubstitution\t0.05\t",
            "chr2\t200\tG\tC\tSubstitution\t0.001\t",
        ]);

        let calls = parse_detection_results(det_file.path()).unwrap();
        let truth = parse_ground_truth(truth_file.path()).unwrap();
        let thresholds = vec![0.0, 0.01, 0.1];

        let report =
            benchmark::run_benchmark(&calls, &truth, &[0.0, 1.0], Some(&thresholds));

        let sweep = report.threshold_sweep.as_ref().unwrap();
        assert_eq!(sweep.len(), 3);
        // At threshold 0.0: both detected
        assert_eq!(sweep[0].sensitivity, 1.0);
        // At threshold 0.01: only chr1 (rvaf=0.05) passes
        assert_eq!(sweep[1].sensitivity, 0.5);
    }

    #[test]
    fn test_chr_normalization_in_benchmark() {
        // Detection uses "1", truth uses "chr1"
        let det_file = write_detection_tsv(&[
            "s1\ttarget1\tSubstitution\tA>T\t0.05\t1.5\t10\t20\tACGT\tTCGT\tvs_ref\t1\t100\tA\tT",
        ]);
        let truth_file = write_truth_tsv(&[
            "chr1\t100\tA\tT\tSubstitution\t0.05\t",
        ]);

        let calls = parse_detection_results(det_file.path()).unwrap();
        let truth = parse_ground_truth(truth_file.path()).unwrap();

        let report = benchmark::run_benchmark(&calls, &truth, &[0.0, 1.0], None);
        assert_eq!(report.confusion.tp, 1);
        assert_eq!(report.confusion.fp, 0);
    }

    #[test]
    fn test_nonexistent_results_file() {
        let result = parse_detection_results(Path::new("/nonexistent/file.tsv"));
        assert!(result.is_err());
    }

    #[test]
    fn test_nonexistent_truth_file() {
        let result = parse_ground_truth(Path::new("/nonexistent/truth.tsv"));
        assert!(result.is_err());
    }
}
