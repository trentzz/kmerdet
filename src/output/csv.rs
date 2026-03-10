use std::io::Write;

use anyhow::Result;

use crate::variant::VariantCall;

pub fn write(calls: &[VariantCall], writer: &mut dyn Write, no_header: bool) -> Result<()> {
    let mut csv_writer = ::csv::WriterBuilder::new().from_writer(writer);

    if !no_header {
        csv_writer.write_record(super::DETECT_HEADERS)?;
    }

    for call in calls {
        csv_writer.write_record(&[
            &call.sample,
            &call.target,
            &call.variant_type.to_string(),
            &call.variant_name,
            &format!("{:.6}", call.rvaf),
            &format!("{:.2}", call.expression),
            &call.min_coverage.to_string(),
            &call.start_kmer_count.to_string(),
            &call.ref_sequence,
            &call.alt_sequence,
            &call.info,
            call.chrom.as_deref().unwrap_or(""),
            &call.pos.map(|p| p.to_string()).unwrap_or_default(),
            call.ref_allele.as_deref().unwrap_or(""),
            call.alt_allele.as_deref().unwrap_or(""),
        ])?;
    }

    csv_writer.flush()?;
    Ok(())
}
