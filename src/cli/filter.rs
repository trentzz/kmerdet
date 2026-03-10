use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

use crate::filter::{ExpectedVariant, FilterResult};
use crate::variant::VariantCall;

#[derive(clap::Args, Debug)]
pub struct FilterArgs {
    /// Detection results input file (from detect/merge)
    #[arg(short, long)]
    pub input: PathBuf,

    /// Expected variants file (TSV: CHROM, POS, REF, ALT, TYPE)
    #[arg(long)]
    pub targets: PathBuf,

    /// Match by full ALT_SEQUENCE instead of POS/REF/ALT
    #[arg(long)]
    pub use_alt: bool,

    /// Min k-mer coverage
    #[arg(long, default_value = "3")]
    pub min_coverage: u32,

    /// Min variant allele frequency
    #[arg(long, default_value = "0.0")]
    pub min_vaf: f64,

    /// Min expression value
    #[arg(long, default_value = "0.0")]
    pub min_expression: f64,

    /// Variant types to include [default: all]
    #[arg(long, num_args = 1..)]
    pub types: Vec<String>,
}

pub fn run(args: FilterArgs, global: &super::GlobalOptions) -> Result<()> {
    // 1. Parse detection results from input TSV file
    let calls = parse_detection_results(&args.input)?;

    // 2. Parse expected variants from --targets TSV file
    let expected = parse_expected_variants(&args.targets)?;

    // 3. Build FilterConfig from CLI args
    let config = crate::filter::FilterConfig {
        min_coverage: args.min_coverage,
        min_vaf: args.min_vaf,
        min_expression: args.min_expression,
        use_alt: args.use_alt,
        types: args.types.clone(),
    };

    // 4. Run filtering
    let results = crate::filter::filter_results(&calls, &expected, &config)?;

    // 5. Write results as TSV
    let mut output: Box<dyn Write> = match &global.output {
        Some(path) => Box::new(
            File::create(path)
                .with_context(|| format!("creating output file: {}", path.display()))?,
        ),
        None => Box::new(std::io::stdout().lock()),
    };

    write_filter_results(&results, &mut output, global.no_header)?;

    Ok(())
}

/// Parse detection results from a TSV file into VariantCall structs.
fn parse_detection_results(path: &Path) -> Result<Vec<VariantCall>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("opening detection results: {}", path.display()))?;

    let mut calls = Vec::new();
    for record in reader.records() {
        let record = record?;
        let call = VariantCall {
            sample: record.get(0).unwrap_or("").to_string(),
            target: record.get(1).unwrap_or("").to_string(),
            variant_type: record.get(2).unwrap_or("Reference").parse()?,
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
        };
        calls.push(call);
    }
    Ok(calls)
}

/// Parse expected variants from a TSV file.
///
/// Expected format: CHROM\tPOS\tREF\tALT\tTYPE (tab-separated, with header)
fn parse_expected_variants(path: &Path) -> Result<Vec<ExpectedVariant>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("opening expected variants: {}", path.display()))?;

    let mut variants = Vec::new();
    for record in reader.records() {
        let record = record?;
        variants.push(ExpectedVariant {
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
    results: &[FilterResult],
    writer: &mut dyn Write,
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

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_write_filter_results_with_header() {
        let results = vec![FilterResult {
            sample: "sample1".to_string(),
            chrom: "chr1".to_string(),
            pos: 12345,
            ref_allele: "A".to_string(),
            alt_allele: "T".to_string(),
            variant_type: "Substitution".to_string(),
            found: "Found".to_string(),
            filter_notes: "PASS".to_string(),
            kmer_vaf: Some(0.05),
            kmer_min_coverage: Some(10),
            kmer_expression: Some(1.5),
            ref_sequence: Some("ACGT".to_string()),
            variant_sequence: Some("TCGT".to_string()),
        }];

        let mut buf = Vec::new();
        write_filter_results(&results, &mut buf, false).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.starts_with("sample\tchrom\tpos\t"));
        assert!(output.contains("sample1\tchr1\t12345\t"));
    }

    #[test]
    fn test_write_filter_results_no_header() {
        let results = vec![FilterResult {
            sample: "s1".to_string(),
            chrom: "chr2".to_string(),
            pos: 100,
            ref_allele: "G".to_string(),
            alt_allele: "C".to_string(),
            variant_type: "Substitution".to_string(),
            found: "Not Found".to_string(),
            filter_notes: "coverage 1 < 3".to_string(),
            kmer_vaf: None,
            kmer_min_coverage: None,
            kmer_expression: None,
            ref_sequence: None,
            variant_sequence: None,
        }];

        let mut buf = Vec::new();
        write_filter_results(&results, &mut buf, true).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(!output.starts_with("sample\t"));
        assert!(output.starts_with("s1\t"));
    }
}
