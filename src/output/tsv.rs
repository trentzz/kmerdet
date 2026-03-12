use std::io::Write;

use anyhow::Result;

use crate::variant::VariantCall;

pub fn write(calls: &[VariantCall], writer: &mut dyn Write, no_header: bool) -> Result<()> {
    if !no_header {
        writeln!(writer, "{}", super::DETECT_HEADERS.join("\t"))?;
    }

    for call in calls {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{:.6}\t{:.2}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            call.sample,
            call.target,
            call.variant_type,
            call.variant_name,
            call.rvaf,
            call.expression,
            call.min_coverage,
            call.start_kmer_count,
            call.ref_sequence,
            call.alt_sequence,
            call.info,
            call.chrom.as_deref().unwrap_or(""),
            call.pos.map(|p| p.to_string()).unwrap_or_default(),
            call.ref_allele.as_deref().unwrap_or(""),
            call.alt_allele.as_deref().unwrap_or(""),
            call.pvalue.map(|p| format!("{:.6e}", p)).unwrap_or_default(),
            call.qual.map(|q| format!("{:.2}", q)).unwrap_or_default(),
            call.ci_lower.map(|v| format!("{:.6}", v)).unwrap_or_default(),
            call.ci_upper.map(|v| format!("{:.6}", v)).unwrap_or_default(),
        )?;
    }

    Ok(())
}
