/// INDEL normalization and deduplication.
///
/// When the k-mer walker discovers deletions or insertions in repeat regions
/// (e.g., deletion of one T from TTTTT), the graph may produce multiple paths
/// representing the same INDEL at different positions. This module provides:
///
/// 1. **Normalization**: Left-align INDELs to a canonical position following
///    VCF normalization conventions (shift left until the alleles differ).
/// 2. **Deduplication**: Collapse variant calls with identical normalized
///    representations, keeping the best-supported call.

use std::collections::HashMap;

use super::{VariantCall, VariantType};

/// A normalized representation of an INDEL variant.
///
/// Two INDELs that represent the same biological event will have identical
/// `NormalizedVariant` values after left-alignment.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct NormalizedVariant {
    /// Chromosome (if available).
    pub chrom: Option<String>,
    /// Left-aligned position within the reference sequence.
    pub pos: usize,
    /// Normalized reference allele.
    pub ref_allele: String,
    /// Normalized alternative allele.
    pub alt_allele: String,
}

/// Normalize an INDEL by left-aligning it within the given reference context.
///
/// Left-alignment shifts the variant as far left as possible while keeping
/// the representation equivalent. This follows VCF normalization conventions:
///
/// 1. While the last bases of ref and alt are identical, trim them.
/// 2. While the first bases of ref and alt are identical, trim and shift left.
/// 3. If both alleles become empty after trimming, the variant was Reference.
///
/// The `ref_context` parameter is the reference DNA sequence that the variant
/// occurs within. `pos` is the 0-based position of the variant's ref_allele
/// within `ref_context`.
///
/// Returns `None` if the variant normalizes to nothing (identical alleles).
pub fn normalize_indel(
    pos: usize,
    ref_allele: &str,
    alt_allele: &str,
    ref_context: &str,
) -> Option<NormalizedVariant> {
    let mut ref_bytes: Vec<u8> = ref_allele.as_bytes().to_vec();
    let mut alt_bytes: Vec<u8> = alt_allele.as_bytes().to_vec();
    let mut current_pos = pos;
    let context_bytes = ref_context.as_bytes();

    // Step 1: Right-trim — remove matching suffix bases.
    while !ref_bytes.is_empty()
        && !alt_bytes.is_empty()
        && ref_bytes.last() == alt_bytes.last()
    {
        ref_bytes.pop();
        alt_bytes.pop();
    }

    // Step 2: Left-trim — remove matching prefix bases, shifting position right.
    while !ref_bytes.is_empty()
        && !alt_bytes.is_empty()
        && ref_bytes[0] == alt_bytes[0]
    {
        ref_bytes.remove(0);
        alt_bytes.remove(0);
        current_pos += 1;
    }

    // If both are empty after trimming, this is not a real variant.
    if ref_bytes.is_empty() && alt_bytes.is_empty() {
        return None;
    }

    // Step 3: Left-align — shift the INDEL left through homopolymer/repeat regions.
    // For a pure insertion (empty ref) or pure deletion (empty alt), shift left
    // while the last base of the non-empty allele matches the base before current_pos
    // in the reference context.
    if ref_bytes.is_empty() {
        // Pure insertion: shift the inserted sequence left.
        while current_pos > 0 {
            let prev_base = context_bytes[current_pos - 1];
            let last_ins = alt_bytes[alt_bytes.len() - 1];
            if prev_base == last_ins {
                // Rotate the inserted sequence: move last base to front.
                let b = alt_bytes.pop().unwrap();
                alt_bytes.insert(0, b);
                current_pos -= 1;
            } else {
                break;
            }
        }
    } else if alt_bytes.is_empty() {
        // Pure deletion: shift the deleted sequence left.
        while current_pos > 0 {
            let prev_base = context_bytes[current_pos - 1];
            let last_del = ref_bytes[ref_bytes.len() - 1];
            if prev_base == last_del {
                // Rotate the deleted sequence: move last base to front.
                let b = ref_bytes.pop().unwrap();
                ref_bytes.insert(0, b);
                current_pos -= 1;
            } else {
                break;
            }
        }
    }
    // For complex INDELs (both non-empty after trimming), no further left-alignment
    // is performed because the alleles are fully trimmed and distinct.

    Some(NormalizedVariant {
        chrom: None,
        pos: current_pos,
        ref_allele: String::from_utf8(ref_bytes).unwrap_or_default(),
        alt_allele: String::from_utf8(alt_bytes).unwrap_or_default(),
    })
}

/// Deduplicate variant calls by collapsing INDELs with equivalent normalized representations.
///
/// SNVs and other non-INDEL types are kept as-is (passed through unchanged).
/// Among duplicate INDELs, the call with the highest expression (best supported)
/// is kept.
///
/// The deduplication key combines (sample, target, normalized position, normalized
/// ref allele, normalized alt allele) so that only truly equivalent variants from
/// the same sample and target are collapsed.
pub fn deduplicate_calls(calls: Vec<VariantCall>) -> Vec<VariantCall> {
    // Separate INDEL calls from non-INDEL calls.
    let mut non_indels: Vec<VariantCall> = Vec::new();
    // Map from dedup key to the best (highest expression) call.
    let mut indel_map: HashMap<IndelDedupKey, VariantCall> = HashMap::new();
    // Track insertion order for deterministic output.
    let mut indel_order: Vec<IndelDedupKey> = Vec::new();

    for call in calls {
        if is_indel_type(call.variant_type) {
            let normalized = normalize_call(&call);
            let key = IndelDedupKey {
                sample: call.sample.clone(),
                target: call.target.clone(),
                normalized,
            };
            match indel_map.get(&key) {
                Some(existing) if existing.expression >= call.expression => {
                    // Existing call has higher or equal expression; skip this one.
                }
                _ => {
                    if !indel_map.contains_key(&key) {
                        indel_order.push(key.clone());
                    }
                    indel_map.insert(key, call);
                }
            }
        } else {
            non_indels.push(call);
        }
    }

    // Reconstruct the result: non-INDELs first, then deduplicated INDELs in order.
    let mut result = non_indels;
    for key in indel_order {
        if let Some(call) = indel_map.remove(&key) {
            result.push(call);
        }
    }

    result
}

/// Check whether a variant type is an INDEL that should be subject to deduplication.
fn is_indel_type(vt: VariantType) -> bool {
    matches!(vt, VariantType::Insertion | VariantType::Deletion | VariantType::Indel)
}

/// Deduplication key for INDEL calls.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct IndelDedupKey {
    sample: String,
    target: String,
    normalized: NormalizedVariant,
}

/// Normalize a variant call for deduplication.
///
/// Uses the call's ref_sequence as the reference context for left-alignment.
/// Falls back to a trivial normalization if position cannot be determined.
fn normalize_call(call: &VariantCall) -> NormalizedVariant {
    let ref_allele = call.ref_allele.as_deref().unwrap_or("");
    let alt_allele = call.alt_allele.as_deref().unwrap_or("");

    // Parse the position from the variant_name ("start:ref/alt:end").
    let pos = parse_start_position(&call.variant_name).unwrap_or(0);

    // Use the ref_sequence as context for left-alignment.
    let ref_context = &call.ref_sequence;

    if ref_context.is_empty() || (ref_allele.is_empty() && alt_allele.is_empty()) {
        // No context or no alleles: return a trivial normalization.
        return NormalizedVariant {
            chrom: call.chrom.clone(),
            pos,
            ref_allele: ref_allele.to_string(),
            alt_allele: alt_allele.to_string(),
        };
    }

    match normalize_indel(pos, ref_allele, alt_allele, ref_context) {
        Some(mut nv) => {
            nv.chrom = call.chrom.clone();
            nv
        }
        None => {
            // Normalization collapsed to nothing; shouldn't happen for real INDELs.
            NormalizedVariant {
                chrom: call.chrom.clone(),
                pos,
                ref_allele: ref_allele.to_string(),
                alt_allele: alt_allele.to_string(),
            }
        }
    }
}

/// Parse the start position from a variant_name string.
///
/// The format is "start:ref/alt:end".
fn parse_start_position(variant_name: &str) -> Option<usize> {
    let first_part = variant_name.split(':').next()?;
    first_part.parse().ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // normalize_indel tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_normalize_deletion_in_homopolymer() {
        // AAAAA: deleting the A at position 3 should left-align to position 0.
        let result = normalize_indel(3, "A", "", "AAAAA").unwrap();
        assert_eq!(result.pos, 0);
        assert_eq!(result.ref_allele, "A");
        assert_eq!(result.alt_allele, "");
    }

    #[test]
    fn test_normalize_deletion_at_different_positions_same_result() {
        // Deleting one T from TTTTT at position 1 vs position 3 should produce
        // the same normalized representation.
        let norm1 = normalize_indel(1, "T", "", "TTTTT").unwrap();
        let norm2 = normalize_indel(3, "T", "", "TTTTT").unwrap();
        assert_eq!(norm1, norm2);
    }

    #[test]
    fn test_normalize_insertion_in_homopolymer() {
        // Inserting an A at position 3 into AAAAA should left-align to position 0.
        let result = normalize_indel(3, "", "A", "AAAAA").unwrap();
        assert_eq!(result.pos, 0);
        assert_eq!(result.ref_allele, "");
        assert_eq!(result.alt_allele, "A");
    }

    #[test]
    fn test_normalize_deletion_no_repeat() {
        // ACGT: deleting C at position 1 should stay at position 1 (no repeat).
        let result = normalize_indel(1, "C", "", "ACGT").unwrap();
        assert_eq!(result.pos, 1);
        assert_eq!(result.ref_allele, "C");
        assert_eq!(result.alt_allele, "");
    }

    #[test]
    fn test_normalize_insertion_no_repeat() {
        // ACGT: inserting G at position 2 should stay at position 2.
        let result = normalize_indel(2, "", "G", "ACGT").unwrap();
        assert_eq!(result.pos, 2);
        assert_eq!(result.ref_allele, "");
        assert_eq!(result.alt_allele, "G");
    }

    #[test]
    fn test_normalize_identical_alleles_returns_none() {
        let result = normalize_indel(2, "AC", "AC", "AACGT");
        assert!(result.is_none());
    }

    #[test]
    fn test_normalize_dinucleotide_repeat() {
        // ATATAT: deleting AT at position 2 vs position 4 should normalize to same position.
        let norm1 = normalize_indel(2, "AT", "", "ATATAT").unwrap();
        let norm2 = normalize_indel(4, "AT", "", "ATATAT").unwrap();
        assert_eq!(norm1, norm2);
    }

    #[test]
    fn test_normalize_complex_indel() {
        // Complex INDEL: both ref and alt non-empty and different.
        // "ACG" -> "TT" at position 2 in context "XXACGXX"
        let result = normalize_indel(2, "ACG", "TT", "XXACGXX").unwrap();
        assert_eq!(result.pos, 2);
        assert_eq!(result.ref_allele, "ACG");
        assert_eq!(result.alt_allele, "TT");
    }

    // -----------------------------------------------------------------------
    // deduplicate_calls tests
    // -----------------------------------------------------------------------

    /// Helper to make a VariantCall for testing.
    fn make_call(
        sample: &str,
        target: &str,
        variant_type: VariantType,
        variant_name: &str,
        ref_allele: &str,
        alt_allele: &str,
        ref_sequence: &str,
        expression: f64,
    ) -> VariantCall {
        VariantCall {
            sample: sample.to_string(),
            target: target.to_string(),
            variant_type,
            variant_name: variant_name.to_string(),
            rvaf: 0.1,
            expression,
            min_coverage: 10,
            start_kmer_count: 10,
            ref_sequence: ref_sequence.to_string(),
            alt_sequence: String::new(),
            info: "vs_ref".to_string(),
            chrom: None,
            pos: None,
            ref_allele: Some(ref_allele.to_string()),
            alt_allele: Some(alt_allele.to_string()),
        }
    }

    #[test]
    fn test_dedup_identical_indels_different_positions_homopolymer() {
        // Two calls representing deletion of one T from TTTTT at different positions.
        // Both should normalize to position 0 and be collapsed.
        let call1 = make_call(
            "sample1", "target1", VariantType::Deletion,
            "1:T/:1", "T", "", "TTTTT", 100.0,
        );
        let call2 = make_call(
            "sample1", "target1", VariantType::Deletion,
            "3:T/:3", "T", "", "TTTTT", 80.0,
        );

        let result = deduplicate_calls(vec![call1, call2]);
        // Should keep only one call (the one with expression 100.0).
        let deletions: Vec<&VariantCall> = result
            .iter()
            .filter(|c| c.variant_type == VariantType::Deletion)
            .collect();
        assert_eq!(deletions.len(), 1);
        assert_eq!(deletions[0].expression, 100.0);
    }

    #[test]
    fn test_dedup_snvs_not_affected() {
        // SNVs should pass through unchanged, even if they look similar.
        let snv1 = make_call(
            "s", "t", VariantType::Substitution,
            "5:A/T:5", "A", "T", "ACGTACGT", 50.0,
        );
        let snv2 = make_call(
            "s", "t", VariantType::Substitution,
            "5:A/T:5", "A", "T", "ACGTACGT", 30.0,
        );

        let result = deduplicate_calls(vec![snv1, snv2]);
        let subs: Vec<&VariantCall> = result
            .iter()
            .filter(|c| c.variant_type == VariantType::Substitution)
            .collect();
        assert_eq!(subs.len(), 2, "SNVs should not be deduplicated");
    }

    #[test]
    fn test_dedup_keeps_highest_expression() {
        // Among duplicates, the highest expression call should be kept.
        let call_low = make_call(
            "s", "t", VariantType::Deletion,
            "2:T/:2", "T", "", "TTTTT", 10.0,
        );
        let call_high = make_call(
            "s", "t", VariantType::Deletion,
            "4:T/:4", "T", "", "TTTTT", 200.0,
        );
        let call_mid = make_call(
            "s", "t", VariantType::Deletion,
            "1:T/:1", "T", "", "TTTTT", 50.0,
        );

        let result = deduplicate_calls(vec![call_low, call_high, call_mid]);
        let deletions: Vec<&VariantCall> = result
            .iter()
            .filter(|c| c.variant_type == VariantType::Deletion)
            .collect();
        assert_eq!(deletions.len(), 1);
        assert_eq!(deletions[0].expression, 200.0);
    }

    #[test]
    fn test_dedup_mixed_indel_and_non_indel() {
        // Mix of types: SNV, Reference, and two duplicate Deletions.
        let ref_call = make_call(
            "s", "t", VariantType::Reference,
            "", "", "", "TTTTT", 1000.0,
        );
        let snv_call = make_call(
            "s", "t", VariantType::Substitution,
            "3:A/T:3", "A", "T", "TTTTT", 50.0,
        );
        let del1 = make_call(
            "s", "t", VariantType::Deletion,
            "1:T/:1", "T", "", "TTTTT", 100.0,
        );
        let del2 = make_call(
            "s", "t", VariantType::Deletion,
            "3:T/:3", "T", "", "TTTTT", 80.0,
        );

        let result = deduplicate_calls(vec![ref_call, snv_call, del1, del2]);

        // Reference and SNV pass through, two deletions become one.
        assert_eq!(result.len(), 3);

        let refs: Vec<_> = result.iter().filter(|c| c.variant_type == VariantType::Reference).collect();
        let subs: Vec<_> = result.iter().filter(|c| c.variant_type == VariantType::Substitution).collect();
        let dels: Vec<_> = result.iter().filter(|c| c.variant_type == VariantType::Deletion).collect();

        assert_eq!(refs.len(), 1);
        assert_eq!(subs.len(), 1);
        assert_eq!(dels.len(), 1);
        assert_eq!(dels[0].expression, 100.0);
    }

    #[test]
    fn test_dedup_no_duplicates_passes_through() {
        // All distinct calls should pass through unchanged.
        let del_t = make_call(
            "s", "t", VariantType::Deletion,
            "2:T/:2", "T", "", "ACGTACGT", 100.0,
        );
        let del_c = make_call(
            "s", "t", VariantType::Deletion,
            "3:G/:3", "G", "", "ACGTACGT", 80.0,
        );
        let ins_a = make_call(
            "s", "t", VariantType::Insertion,
            "4:/A:4", "", "A", "ACGTACGT", 60.0,
        );

        let result = deduplicate_calls(vec![del_t, del_c, ins_a]);
        assert_eq!(result.len(), 3, "no duplicates should mean no collapsing");
    }

    #[test]
    fn test_dedup_different_samples_not_collapsed() {
        // Same INDEL in different samples should NOT be collapsed.
        let call1 = make_call(
            "sample_A", "t", VariantType::Deletion,
            "1:T/:1", "T", "", "TTTTT", 100.0,
        );
        let call2 = make_call(
            "sample_B", "t", VariantType::Deletion,
            "3:T/:3", "T", "", "TTTTT", 80.0,
        );

        let result = deduplicate_calls(vec![call1, call2]);
        assert_eq!(result.len(), 2, "different samples should not be collapsed");
    }

    #[test]
    fn test_dedup_different_targets_not_collapsed() {
        // Same INDEL in different targets should NOT be collapsed.
        let call1 = make_call(
            "s", "target_A", VariantType::Deletion,
            "1:T/:1", "T", "", "TTTTT", 100.0,
        );
        let call2 = make_call(
            "s", "target_B", VariantType::Deletion,
            "3:T/:3", "T", "", "TTTTT", 80.0,
        );

        let result = deduplicate_calls(vec![call1, call2]);
        assert_eq!(result.len(), 2, "different targets should not be collapsed");
    }

    #[test]
    fn test_dedup_empty_input() {
        let result = deduplicate_calls(vec![]);
        assert!(result.is_empty());
    }

    #[test]
    fn test_dedup_insertion_duplicates() {
        // Two insertions of A at different positions in a homopolymer.
        let call1 = make_call(
            "s", "t", VariantType::Insertion,
            "1:/A:1", "", "A", "AAAA", 120.0,
        );
        let call2 = make_call(
            "s", "t", VariantType::Insertion,
            "3:/A:3", "", "A", "AAAA", 90.0,
        );

        let result = deduplicate_calls(vec![call1, call2]);
        let insertions: Vec<&VariantCall> = result
            .iter()
            .filter(|c| c.variant_type == VariantType::Insertion)
            .collect();
        assert_eq!(insertions.len(), 1);
        assert_eq!(insertions[0].expression, 120.0);
    }

    #[test]
    fn test_dedup_complex_indel_duplicates() {
        // Two complex INDELs with the same normalized representation.
        // Position 2: "ACG" -> "TT" — same alleles, same position, should dedup.
        let call1 = make_call(
            "s", "t", VariantType::Indel,
            "2:ACG/TT:4", "ACG", "TT", "XXACGXX", 50.0,
        );
        let call2 = make_call(
            "s", "t", VariantType::Indel,
            "2:ACG/TT:4", "ACG", "TT", "XXACGXX", 70.0,
        );

        let result = deduplicate_calls(vec![call1, call2]);
        let indels: Vec<&VariantCall> = result
            .iter()
            .filter(|c| c.variant_type == VariantType::Indel)
            .collect();
        assert_eq!(indels.len(), 1);
        assert_eq!(indels[0].expression, 70.0);
    }
}
