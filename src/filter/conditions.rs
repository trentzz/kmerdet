/// Filter conditions: type match, coverage, VAF thresholds.
///
/// Each condition is a predicate that a variant call must pass.

use crate::variant::VariantCall;

/// Check if a call meets the minimum coverage threshold.
pub fn passes_coverage(call: &VariantCall, min_coverage: u32) -> bool {
    call.min_coverage >= min_coverage as u64
}

/// Check if a call meets the minimum VAF threshold.
pub fn passes_vaf(call: &VariantCall, min_vaf: f64) -> bool {
    call.rvaf >= min_vaf
}

/// Check if a call meets the minimum expression threshold.
pub fn passes_expression(call: &VariantCall, min_expression: f64) -> bool {
    call.expression >= min_expression
}

/// Check if a call's variant type is in the allowed types list.
/// Empty list means all types are allowed.
pub fn passes_type(call: &VariantCall, types: &[String]) -> bool {
    if types.is_empty() {
        return true;
    }
    let call_type = call.variant_type.to_string().to_lowercase();
    types.iter().any(|t| t.to_lowercase() == call_type)
}
