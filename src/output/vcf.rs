use std::io::Write;

use anyhow::Result;

use crate::variant::VariantCall;

pub fn write(calls: &[VariantCall], writer: &mut dyn Write) -> Result<()> {
    // VCF 4.3 header
    writeln!(writer, "##fileformat=VCFv4.3")?;
    writeln!(
        writer,
        "##source=kmerdet v{}",
        env!("CARGO_PKG_VERSION")
    )?;
    writeln!(
        writer,
        "##INFO=<ID=KVAF,Number=1,Type=Float,Description=\"K-mer variant allele frequency\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=KCOV,Number=1,Type=Integer,Description=\"Minimum k-mer coverage along path\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=KEXP,Number=1,Type=Float,Description=\"K-mer expression estimate\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=KTYPE,Number=1,Type=String,Description=\"Variant type from k-mer analysis\">"
    )?;
    writeln!(
        writer,
        "##INFO=<ID=KPV,Number=1,Type=Float,Description=\"K-mer variant p-value\">"
    )?;
    writeln!(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;

    for call in calls {
        let chrom = call.chrom.as_deref().unwrap_or(".");
        let pos = call.pos.unwrap_or(0);
        let ref_allele = call.ref_allele.as_deref().unwrap_or(".");
        let alt_allele = call.alt_allele.as_deref().unwrap_or(".");
        let qual_str = call
            .qual
            .map(|q| format!("{:.1}", q))
            .unwrap_or_else(|| ".".to_string());

        let mut info = format!(
            "KVAF={:.6};KCOV={};KEXP={:.2};KTYPE={}",
            call.rvaf, call.min_coverage, call.expression, call.variant_type,
        );
        if let Some(pv) = call.pvalue {
            info.push_str(&format!(";KPV={:.6e}", pv));
        }

        writeln!(
            writer,
            "{}\t{}\t.\t{}\t{}\t{}\tPASS\t{}",
            chrom, pos, ref_allele, alt_allele, qual_str, info,
        )?;
    }

    Ok(())
}
