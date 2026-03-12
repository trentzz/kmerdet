pub mod csv;
pub mod excel;
pub mod json;
pub mod tsv;
pub mod vcf;

use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result};

use crate::variant::VariantCall;

/// Write variant calls to the given writer in the specified format.
pub fn write_calls(
    calls: &[VariantCall],
    format: &crate::cli::OutputFormat,
    writer: &mut dyn Write,
    no_header: bool,
) -> Result<()> {
    match format {
        crate::cli::OutputFormat::Tsv => tsv::write(calls, writer, no_header),
        crate::cli::OutputFormat::Csv => csv::write(calls, writer, no_header),
        crate::cli::OutputFormat::Vcf => vcf::write(calls, writer),
        crate::cli::OutputFormat::Json => json::write_json(calls, writer),
        crate::cli::OutputFormat::Jsonl => json::write_jsonl(calls, writer),
        crate::cli::OutputFormat::Xlsx => {
            anyhow::bail!("XLSX output requires a file path, not a stream writer")
        }
    }
}

/// Write variant calls to an XLSX file.
pub fn write_calls_xlsx(
    calls: &[VariantCall],
    path: &std::path::Path,
) -> Result<()> {
    excel::write(calls, path)
}

/// Column headers for detection output.
pub const DETECT_HEADERS: &[&str] = &[
    "sample",
    "target",
    "type",
    "variant_name",
    "rVAF",
    "expression",
    "min_coverage",
    "start_kmer_count",
    "ref_sequence",
    "alt_sequence",
    "info",
    "chrom",
    "pos",
    "ref_allele",
    "alt_allele",
];

/// Parse detection results from a TSV file into VariantCall structs.
///
/// Shared utility used by filter, stats, plot, and benchmark subcommands.
pub fn parse_detection_tsv(path: &Path) -> Result<Vec<VariantCall>> {
    let mut reader = ::csv::ReaderBuilder::new()
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
            path_score: record.get(6).unwrap_or("0").parse().unwrap_or(0),
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
