use std::io::Write;

use anyhow::Result;

use crate::variant::VariantCall;

/// Write all calls as a single JSON array.
pub fn write_json(calls: &[VariantCall], writer: &mut dyn Write) -> Result<()> {
    serde_json::to_writer_pretty(writer, calls)?;
    Ok(())
}

/// Write calls as newline-delimited JSON (one JSON object per line).
pub fn write_jsonl(calls: &[VariantCall], writer: &mut dyn Write) -> Result<()> {
    for call in calls {
        serde_json::to_writer(&mut *writer, call)?;
        writeln!(writer)?;
    }
    Ok(())
}
