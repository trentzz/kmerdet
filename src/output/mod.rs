pub mod csv;
pub mod excel;
pub mod json;
pub mod tsv;
pub mod vcf;

use std::io::Write;

use anyhow::Result;

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
