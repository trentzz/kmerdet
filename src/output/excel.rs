use std::path::Path;

use anyhow::Result;
use rust_xlsxwriter::{Format, Workbook};

use crate::variant::VariantCall;

pub fn write(calls: &[VariantCall], path: &Path) -> Result<()> {
    let mut workbook = Workbook::new();
    let worksheet = workbook.add_worksheet();

    // Header formatting
    let header_fmt = Format::new().set_bold();

    // Write headers
    for (col, &header) in super::DETECT_HEADERS.iter().enumerate() {
        worksheet.write_string_with_format(0, col as u16, header, &header_fmt)?;
    }

    // Write data rows
    for (row_idx, call) in calls.iter().enumerate() {
        let row = (row_idx + 1) as u32;
        worksheet.write_string(row, 0, &call.sample)?;
        worksheet.write_string(row, 1, &call.target)?;
        worksheet.write_string(row, 2, &call.variant_type.to_string())?;
        worksheet.write_string(row, 3, &call.variant_name)?;
        worksheet.write_number(row, 4, call.rvaf)?;
        worksheet.write_number(row, 5, call.expression)?;
        worksheet.write_number(row, 6, call.min_coverage as f64)?;
        worksheet.write_number(row, 7, call.start_kmer_count as f64)?;
        worksheet.write_string(row, 8, &call.ref_sequence)?;
        worksheet.write_string(row, 9, &call.alt_sequence)?;
        worksheet.write_string(row, 10, &call.info)?;
        worksheet.write_string(row, 11, call.chrom.as_deref().unwrap_or(""))?;
        if let Some(pos) = call.pos {
            worksheet.write_number(row, 12, pos as f64)?;
        }
        worksheet.write_string(row, 13, call.ref_allele.as_deref().unwrap_or(""))?;
        worksheet.write_string(row, 14, call.alt_allele.as_deref().unwrap_or(""))?;
    }

    workbook.save(path)?;
    Ok(())
}
