/// Multi-k consensus voting for variant detection.
///
/// When detection runs independently at multiple k-mer lengths (e.g., k=21, 31, 41),
/// this module merges the results using weighted consensus voting. Shorter k values
/// improve INDEL sensitivity while longer k values maintain SNV specificity.
///
/// Variants from different k values are matched by (sample, target, normalized variant)
/// and merged with coverage-weighted rVAF averaging and confidence tiering.

use std::collections::HashMap;

use super::normalize::{normalize_indel, NormalizedVariant};
use super::{VariantCall, VariantType};

/// Confidence tier based on how many k values detected a variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ConfidenceTier {
    /// Detected at all k values — highest confidence.
    Tier1,
    /// Detected at 2 out of 3 (or majority) k values.
    Tier2,
    /// Detected at k=31 only (standard k, moderate confidence).
    Tier3,
    /// Detected at a single non-standard k only — lowest confidence.
    Tier4,
}

impl std::fmt::Display for ConfidenceTier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Tier1 => write!(f, "Tier1"),
            Self::Tier2 => write!(f, "Tier2"),
            Self::Tier3 => write!(f, "Tier3"),
            Self::Tier4 => write!(f, "Tier4"),
        }
    }
}

/// A merged variant call from multiple k values.
#[derive(Debug, Clone)]
pub struct ConsensusCall {
    /// The best representative VariantCall (highest expression among k values).
    pub call: VariantCall,
    /// Which k values detected this variant.
    pub k_values: Vec<u8>,
    /// Coverage-weighted average rVAF.
    pub consensus_rvaf: f64,
    /// Confidence tier.
    pub tier: ConfidenceTier,
    /// Per-k rVAFs.
    pub per_k_rvafs: HashMap<u8, f64>,
    /// Weighted consensus score (0.0–1.0).
    pub consensus_score: f64,
}

/// Compute the consensus weight for a variant type at a given k value.
///
/// Weights reflect domain knowledge: shorter k values are better for INDELs
/// (fewer k-mers span the breakpoint), while longer k values are better for
/// SNVs (more unique context reduces false positives).
pub fn consensus_weight(variant_type: VariantType, k: u8) -> f64 {
    match variant_type {
        VariantType::Substitution => match k {
            k if k <= 21 => 0.7,
            k if k <= 31 => 1.0,
            _ => 0.9,
        },
        VariantType::Insertion | VariantType::Deletion => match k {
            k if k <= 21 => 1.0,
            k if k <= 31 => 0.9,
            _ => 0.7,
        },
        VariantType::Indel => match k {
            k if k <= 21 => 1.0,
            k if k <= 31 => 0.7,
            _ => 0.4,
        },
        VariantType::Itd => match k {
            k if k <= 25 => 1.0,
            k if k <= 31 => 0.8,
            _ => 0.3,
        },
        VariantType::Reference => 1.0,
    }
}

/// Key for grouping variant calls across k values.
///
/// Two calls from different k values represent the same variant if they have
/// the same sample, target, and either:
/// 1. After INDEL normalization, the same (chrom, pos, ref, alt), OR
/// 2. The same variant_name
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct VariantKey {
    sample: String,
    target: String,
    /// The matching identity: normalized coordinates or variant name.
    identity: VariantIdentity,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
enum VariantIdentity {
    /// Matched by normalized variant coordinates.
    Normalized(NormalizedVariant),
    /// Matched by variant name (for SNVs and cases without allele info).
    ByName(String),
    /// Reference calls are grouped by sample+target only.
    Reference,
}

/// Check if two variant types are compatible for merging across k values.
///
/// INDELs can be typed differently at different k values (e.g., Insertion at k=21
/// might be classified as Indel at k=41 due to different breakpoint resolution).
fn types_compatible(a: VariantType, b: VariantType) -> bool {
    if a == b {
        return true;
    }
    // INDEL-family types are compatible with each other.
    let indel_family = |t: VariantType| {
        matches!(
            t,
            VariantType::Insertion | VariantType::Deletion | VariantType::Indel
        )
    };
    indel_family(a) && indel_family(b)
}

/// Extract a VariantKey for grouping across k values.
fn variant_key(call: &VariantCall) -> VariantKey {
    if call.variant_type == VariantType::Reference {
        return VariantKey {
            sample: call.sample.clone(),
            target: call.target.clone(),
            identity: VariantIdentity::Reference,
        };
    }

    // Try to use normalized coordinates for INDEL-family variants.
    if is_indel_type(call.variant_type) {
        if let (Some(ref_al), Some(alt_al)) = (&call.ref_allele, &call.alt_allele) {
            if !ref_al.is_empty() || !alt_al.is_empty() {
                let pos = parse_start_position(&call.variant_name).unwrap_or(0);
                let ref_context = if call.ref_sequence.is_empty() {
                    None
                } else {
                    Some(call.ref_sequence.as_str())
                };
                let normalized = normalize_indel(
                    call.chrom.as_deref().unwrap_or(""),
                    pos as u64,
                    ref_al,
                    alt_al,
                    ref_context,
                );
                return VariantKey {
                    sample: call.sample.clone(),
                    target: call.target.clone(),
                    identity: VariantIdentity::Normalized(normalized),
                };
            }
        }
    }

    // For SNVs and cases without allele info, use variant_name matching.
    // But first, try normalized coordinates for SNVs too if allele info is available.
    if let (Some(ref_al), Some(alt_al)) = (&call.ref_allele, &call.alt_allele) {
        if !ref_al.is_empty() || !alt_al.is_empty() {
            let pos = parse_start_position(&call.variant_name).unwrap_or(0);
            let normalized = normalize_indel(
                call.chrom.as_deref().unwrap_or(""),
                pos as u64,
                ref_al,
                alt_al,
                None, // No ref context normalization for SNVs
            );
            return VariantKey {
                sample: call.sample.clone(),
                target: call.target.clone(),
                identity: VariantIdentity::Normalized(normalized),
            };
        }
    }

    // Fallback: use variant name.
    VariantKey {
        sample: call.sample.clone(),
        target: call.target.clone(),
        identity: VariantIdentity::ByName(call.variant_name.clone()),
    }
}

fn is_indel_type(vt: VariantType) -> bool {
    matches!(
        vt,
        VariantType::Insertion | VariantType::Deletion | VariantType::Indel
    )
}

/// Parse the start position from a variant_name string (format: "start:ref/alt:end").
fn parse_start_position(variant_name: &str) -> Option<usize> {
    let first_part = variant_name.split(':').next()?;
    first_part.parse().ok()
}

/// Assign a confidence tier based on which k values detected the variant.
fn assign_tier(k_values: &[u8], total_k_values: usize) -> ConfidenceTier {
    let detected = k_values.len();
    if detected >= total_k_values && total_k_values > 0 {
        ConfidenceTier::Tier1
    } else if detected >= 2 {
        ConfidenceTier::Tier2
    } else if detected == 1 && k_values[0] == 31 {
        ConfidenceTier::Tier3
    } else {
        ConfidenceTier::Tier4
    }
}

/// Merge variant calls from multiple k values into consensus calls.
///
/// Groups calls by (sample, target, normalized variant), applies weighted
/// consensus voting, computes coverage-weighted rVAF, and assigns confidence tiers.
///
/// # Arguments
/// * `calls_per_k` - Pairs of (k_value, calls_at_that_k).
/// * `total_k_values` - Total number of k values that were run (for tier assignment).
///
/// # Returns
/// Consensus calls sorted by sample, then target.
pub fn merge_multi_k(
    calls_per_k: &[(u8, Vec<VariantCall>)],
    total_k_values: usize,
) -> Vec<ConsensusCall> {
    if calls_per_k.is_empty() {
        return Vec::new();
    }

    // Group calls by VariantKey. For each group, track (k, call) pairs.
    // We need stable ordering, so use a Vec of keys in insertion order.
    let mut groups: HashMap<VariantKey, Vec<(u8, VariantCall)>> = HashMap::new();
    let mut key_order: Vec<VariantKey> = Vec::new();

    for (k, calls) in calls_per_k {
        for call in calls {
            let key = variant_key(call);

            // Check type compatibility: only merge if types are compatible
            // with existing entries in this group.
            let should_merge = if let Some(existing) = groups.get(&key) {
                existing
                    .iter()
                    .all(|(_, ec)| types_compatible(ec.variant_type, call.variant_type))
            } else {
                true
            };

            if should_merge {
                if !groups.contains_key(&key) {
                    key_order.push(key.clone());
                }
                groups.entry(key).or_default().push((*k, call.clone()));
            } else {
                // Type incompatibility — treat as a separate variant.
                // This shouldn't happen often since we use normalized coordinates.
                let fallback_key = VariantKey {
                    sample: call.sample.clone(),
                    target: call.target.clone(),
                    identity: VariantIdentity::ByName(format!(
                        "k{}:{}",
                        k, call.variant_name
                    )),
                };
                if !groups.contains_key(&fallback_key) {
                    key_order.push(fallback_key.clone());
                }
                groups
                    .entry(fallback_key)
                    .or_default()
                    .push((*k, call.clone()));
            }
        }
    }

    // Build ConsensusCall for each group.
    let mut consensus_calls: Vec<ConsensusCall> = Vec::new();

    for key in &key_order {
        let entries = match groups.get(key) {
            Some(e) => e,
            None => continue,
        };

        if entries.is_empty() {
            continue;
        }

        // Determine k values that detected this variant.
        let mut k_values: Vec<u8> = entries.iter().map(|(k, _)| *k).collect();
        k_values.sort_unstable();
        k_values.dedup();

        // Per-k rVAFs.
        let mut per_k_rvafs: HashMap<u8, f64> = HashMap::new();
        for (k, call) in entries {
            // If multiple calls at the same k (shouldn't happen normally),
            // keep the one with higher rVAF.
            let entry = per_k_rvafs.entry(*k).or_insert(0.0);
            if call.rvaf > *entry {
                *entry = call.rvaf;
            }
        }

        // Use the variant type from the best (highest expression) call
        // to compute consensus weights.
        let best_call = entries
            .iter()
            .max_by(|(_, a), (_, b)| {
                a.expression
                    .partial_cmp(&b.expression)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(_, c)| c)
            .unwrap();

        let variant_type = best_call.variant_type;

        // Compute coverage-weighted average rVAF.
        let mut weighted_rvaf_sum = 0.0;
        let mut weight_sum = 0.0;
        for (k, call) in entries {
            let w = consensus_weight(variant_type, *k) * call.min_coverage as f64;
            weighted_rvaf_sum += w * call.rvaf;
            weight_sum += w;
        }
        let consensus_rvaf = if weight_sum > 0.0 {
            weighted_rvaf_sum / weight_sum
        } else {
            best_call.rvaf
        };

        // Compute consensus score: sum of weights for detected k values
        // divided by max possible weight sum.
        let all_k_values: Vec<u8> = calls_per_k.iter().map(|(k, _)| *k).collect();
        let max_weight: f64 = all_k_values
            .iter()
            .map(|k| consensus_weight(variant_type, *k))
            .sum();
        let actual_weight: f64 = k_values
            .iter()
            .map(|k| consensus_weight(variant_type, *k))
            .sum();
        let consensus_score = if max_weight > 0.0 {
            actual_weight / max_weight
        } else {
            0.0
        };

        // Assign confidence tier.
        let tier = assign_tier(&k_values, total_k_values);

        consensus_calls.push(ConsensusCall {
            call: best_call.clone(),
            k_values,
            consensus_rvaf,
            tier,
            per_k_rvafs,
            consensus_score,
        });
    }

    // Sort by sample, then target.
    consensus_calls.sort_by(|a, b| {
        a.call
            .sample
            .cmp(&b.call.sample)
            .then_with(|| a.call.target.cmp(&b.call.target))
    });

    consensus_calls
}

/// Convert a ConsensusCall back into a VariantCall for output.
///
/// Copies the best representative call and enriches the `info` field with
/// consensus metadata (tier, k values, per-k rVAFs, consensus score).
pub fn consensus_to_variant_call(consensus: &ConsensusCall) -> VariantCall {
    let mut call = consensus.call.clone();

    // Update rVAF to the consensus-weighted value.
    call.rvaf = consensus.consensus_rvaf;

    // Enrich the info field with consensus metadata.
    let k_str: Vec<String> = consensus.k_values.iter().map(|k| k.to_string()).collect();
    let per_k_str: Vec<String> = consensus
        .per_k_rvafs
        .iter()
        .map(|(k, v)| format!("k{}={:.4}", k, v))
        .collect();

    call.info = format!(
        "multi_k={};tier={};score={:.3};per_k=[{}]",
        k_str.join(","),
        consensus.tier,
        consensus.consensus_score,
        per_k_str.join(","),
    );

    call
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper to create a VariantCall for testing.
    fn make_call(
        sample: &str,
        target: &str,
        variant_type: VariantType,
        variant_name: &str,
        rvaf: f64,
        expression: f64,
        min_coverage: u64,
        ref_allele: &str,
        alt_allele: &str,
        ref_sequence: &str,
    ) -> VariantCall {
        VariantCall {
            sample: sample.to_string(),
            target: target.to_string(),
            variant_type,
            variant_name: variant_name.to_string(),
            rvaf,
            expression,
            min_coverage,
            start_kmer_count: min_coverage,
            ref_sequence: ref_sequence.to_string(),
            alt_sequence: String::new(),
            info: "vs_ref".to_string(),
            chrom: None,
            pos: None,
            ref_allele: Some(ref_allele.to_string()),
            alt_allele: Some(alt_allele.to_string()),
            pvalue: None,
            qual: None,
            ci_lower: None,
            ci_upper: None,
        }
    }

    // -- consensus_weight tests --

    #[test]
    fn test_consensus_weight_snv_prefers_k31() {
        assert_eq!(consensus_weight(VariantType::Substitution, 21), 0.7);
        assert_eq!(consensus_weight(VariantType::Substitution, 31), 1.0);
        assert_eq!(consensus_weight(VariantType::Substitution, 41), 0.9);
    }

    #[test]
    fn test_consensus_weight_insertion_prefers_short_k() {
        assert_eq!(consensus_weight(VariantType::Insertion, 21), 1.0);
        assert_eq!(consensus_weight(VariantType::Insertion, 31), 0.9);
        assert_eq!(consensus_weight(VariantType::Insertion, 41), 0.7);
    }

    #[test]
    fn test_consensus_weight_deletion_prefers_short_k() {
        assert_eq!(consensus_weight(VariantType::Deletion, 21), 1.0);
        assert_eq!(consensus_weight(VariantType::Deletion, 31), 0.9);
        assert_eq!(consensus_weight(VariantType::Deletion, 41), 0.7);
    }

    #[test]
    fn test_consensus_weight_indel_strongly_prefers_short_k() {
        assert_eq!(consensus_weight(VariantType::Indel, 21), 1.0);
        assert_eq!(consensus_weight(VariantType::Indel, 31), 0.7);
        assert_eq!(consensus_weight(VariantType::Indel, 41), 0.4);
    }

    #[test]
    fn test_consensus_weight_itd() {
        assert_eq!(consensus_weight(VariantType::Itd, 21), 1.0);
        assert_eq!(consensus_weight(VariantType::Itd, 25), 1.0);
        assert_eq!(consensus_weight(VariantType::Itd, 31), 0.8);
        assert_eq!(consensus_weight(VariantType::Itd, 41), 0.3);
    }

    #[test]
    fn test_consensus_weight_reference_always_1() {
        assert_eq!(consensus_weight(VariantType::Reference, 21), 1.0);
        assert_eq!(consensus_weight(VariantType::Reference, 31), 1.0);
        assert_eq!(consensus_weight(VariantType::Reference, 41), 1.0);
    }

    // -- merge_multi_k: Tier1 (all k values detect the same variant) --

    #[test]
    fn test_merge_all_k_values_tier1() {
        let call_k21 = make_call(
            "sample1", "target1", VariantType::Substitution,
            "5:A/T:5", 0.10, 100.0, 50,
            "A", "T", "ACGTACGT",
        );
        let call_k31 = make_call(
            "sample1", "target1", VariantType::Substitution,
            "5:A/T:5", 0.12, 120.0, 60,
            "A", "T", "ACGTACGT",
        );
        let call_k41 = make_call(
            "sample1", "target1", VariantType::Substitution,
            "5:A/T:5", 0.08, 80.0, 40,
            "A", "T", "ACGTACGT",
        );

        let calls_per_k = vec![
            (21, vec![call_k21]),
            (31, vec![call_k31]),
            (41, vec![call_k41]),
        ];

        let result = merge_multi_k(&calls_per_k, 3);
        assert_eq!(result.len(), 1);

        let consensus = &result[0];
        assert_eq!(consensus.tier, ConfidenceTier::Tier1);
        assert_eq!(consensus.k_values, vec![21, 31, 41]);
        assert_eq!(consensus.per_k_rvafs.len(), 3);
        assert!((consensus.per_k_rvafs[&21] - 0.10).abs() < 1e-6);
        assert!((consensus.per_k_rvafs[&31] - 0.12).abs() < 1e-6);
        assert!((consensus.per_k_rvafs[&41] - 0.08).abs() < 1e-6);
        // Best call should be k=31 (highest expression: 120)
        assert_eq!(consensus.call.expression, 120.0);
        // Consensus score should be 1.0 (all k values detected)
        assert!((consensus.consensus_score - 1.0).abs() < 1e-6);
    }

    // -- merge_multi_k: Tier2 (2 out of 3 k values) --

    #[test]
    fn test_merge_two_of_three_tier2() {
        let call_k21 = make_call(
            "sample1", "target1", VariantType::Deletion,
            "5:AT/A:6", 0.15, 90.0, 45,
            "AT", "A", "ACGTATCGT",
        );
        let call_k31 = make_call(
            "sample1", "target1", VariantType::Deletion,
            "5:AT/A:6", 0.18, 110.0, 55,
            "AT", "A", "ACGTATCGT",
        );

        let calls_per_k = vec![
            (21, vec![call_k21]),
            (31, vec![call_k31]),
            (41, vec![]), // k=41 didn't detect it
        ];

        let result = merge_multi_k(&calls_per_k, 3);
        assert_eq!(result.len(), 1);

        let consensus = &result[0];
        assert_eq!(consensus.tier, ConfidenceTier::Tier2);
        assert_eq!(consensus.k_values, vec![21, 31]);
    }

    // -- merge_multi_k: Tier3 (k=31 only) --

    #[test]
    fn test_merge_k31_only_tier3() {
        let call_k31 = make_call(
            "sample1", "target1", VariantType::Substitution,
            "5:A/T:5", 0.05, 50.0, 25,
            "A", "T", "ACGTACGT",
        );

        let calls_per_k = vec![
            (21, vec![]),
            (31, vec![call_k31]),
            (41, vec![]),
        ];

        let result = merge_multi_k(&calls_per_k, 3);
        assert_eq!(result.len(), 1);

        let consensus = &result[0];
        assert_eq!(consensus.tier, ConfidenceTier::Tier3);
        assert_eq!(consensus.k_values, vec![31]);
    }

    // -- merge_multi_k: Tier4 (single non-standard k) --

    #[test]
    fn test_merge_single_nonstandard_k_tier4() {
        let call_k21 = make_call(
            "sample1", "target1", VariantType::Insertion,
            "5:A/AT:5", 0.03, 30.0, 15,
            "A", "AT", "ACGTACGT",
        );

        let calls_per_k = vec![
            (21, vec![call_k21]),
            (31, vec![]),
            (41, vec![]),
        ];

        let result = merge_multi_k(&calls_per_k, 3);
        assert_eq!(result.len(), 1);

        let consensus = &result[0];
        assert_eq!(consensus.tier, ConfidenceTier::Tier4);
        assert_eq!(consensus.k_values, vec![21]);
    }

    // -- Coverage-weighted rVAF averaging --

    #[test]
    fn test_coverage_weighted_rvaf() {
        // k=21: rVAF=0.10, min_coverage=100, weight(SNV,21)=0.7
        // k=31: rVAF=0.20, min_coverage=200, weight(SNV,31)=1.0
        // weighted_rvaf = (0.7*100*0.10 + 1.0*200*0.20) / (0.7*100 + 1.0*200)
        //               = (7 + 40) / (70 + 200) = 47/270 ≈ 0.1741
        let call_k21 = make_call(
            "s", "t", VariantType::Substitution,
            "5:A/T:5", 0.10, 100.0, 100,
            "A", "T", "ACGTACGT",
        );
        let call_k31 = make_call(
            "s", "t", VariantType::Substitution,
            "5:A/T:5", 0.20, 200.0, 200,
            "A", "T", "ACGTACGT",
        );

        let calls_per_k = vec![
            (21, vec![call_k21]),
            (31, vec![call_k31]),
        ];

        let result = merge_multi_k(&calls_per_k, 2);
        assert_eq!(result.len(), 1);

        let expected_rvaf = 47.0 / 270.0;
        assert!(
            (result[0].consensus_rvaf - expected_rvaf).abs() < 1e-4,
            "expected ~{:.4}, got {:.4}",
            expected_rvaf,
            result[0].consensus_rvaf,
        );
    }

    // -- Different variants at same target should NOT merge --

    #[test]
    fn test_different_variants_same_target_no_merge() {
        let snv = make_call(
            "s", "t", VariantType::Substitution,
            "5:A/T:5", 0.10, 100.0, 50,
            "A", "T", "ACGTACGT",
        );
        let del = make_call(
            "s", "t", VariantType::Deletion,
            "3:GT/G:4", 0.15, 90.0, 45,
            "GT", "G", "ACGTACGT",
        );

        let calls_per_k = vec![
            (31, vec![snv, del]),
        ];

        let result = merge_multi_k(&calls_per_k, 1);
        assert_eq!(
            result.len(),
            2,
            "different variants at same target should not merge"
        );
    }

    // -- INDEL matching across k values (different breakpoints normalized) --

    #[test]
    fn test_indel_matching_across_k_normalized() {
        // Same deletion in a homopolymer, but reported at different positions by
        // different k values. After normalization, they should match.
        // ref_sequence = "TTTTTTTTT" (poly-T)
        // k=21: pos 4, T/TT → normalized to pos 0, T/TT (left-aligned)
        // k=31: pos 6, T/TT → normalized to pos 0, T/TT (left-aligned)
        let call_k21 = make_call(
            "s", "t", VariantType::Insertion,
            "4:T/TT:4", 0.10, 100.0, 50,
            "T", "TT", "TTTTTTTTT",
        );
        let call_k31 = make_call(
            "s", "t", VariantType::Insertion,
            "6:T/TT:6", 0.12, 120.0, 60,
            "T", "TT", "TTTTTTTTT",
        );

        let calls_per_k = vec![
            (21, vec![call_k21]),
            (31, vec![call_k31]),
        ];

        let result = merge_multi_k(&calls_per_k, 2);
        assert_eq!(
            result.len(),
            1,
            "same INDEL at different positions should merge after normalization"
        );
        assert_eq!(result[0].k_values, vec![21, 31]);
    }

    // -- Empty input --

    #[test]
    fn test_merge_empty_input() {
        let result = merge_multi_k(&[], 0);
        assert!(result.is_empty());
    }

    #[test]
    fn test_merge_all_empty_call_lists() {
        let calls_per_k: Vec<(u8, Vec<VariantCall>)> = vec![
            (21, vec![]),
            (31, vec![]),
            (41, vec![]),
        ];
        let result = merge_multi_k(&calls_per_k, 3);
        assert!(result.is_empty());
    }

    // -- Reference calls merge --

    #[test]
    fn test_reference_calls_merge() {
        let ref_k21 = make_call(
            "s", "t", VariantType::Reference,
            "", 1.0, 500.0, 250,
            "", "", "ACGTACGT",
        );
        let ref_k31 = make_call(
            "s", "t", VariantType::Reference,
            "", 1.0, 600.0, 300,
            "", "", "ACGTACGT",
        );

        let calls_per_k = vec![
            (21, vec![ref_k21]),
            (31, vec![ref_k31]),
        ];

        let result = merge_multi_k(&calls_per_k, 2);
        // Reference calls at same sample+target should merge.
        let ref_calls: Vec<_> = result
            .iter()
            .filter(|c| c.call.variant_type == VariantType::Reference)
            .collect();
        assert_eq!(ref_calls.len(), 1);
    }

    // -- consensus_to_variant_call --

    #[test]
    fn test_consensus_to_variant_call_preserves_fields() {
        let call = make_call(
            "sample1", "target1", VariantType::Substitution,
            "5:A/T:5", 0.10, 100.0, 50,
            "A", "T", "ACGTACGT",
        );
        let consensus = ConsensusCall {
            call: call.clone(),
            k_values: vec![21, 31, 41],
            consensus_rvaf: 0.12,
            tier: ConfidenceTier::Tier1,
            per_k_rvafs: {
                let mut m = HashMap::new();
                m.insert(21, 0.10);
                m.insert(31, 0.12);
                m.insert(41, 0.14);
                m
            },
            consensus_score: 1.0,
        };

        let result = consensus_to_variant_call(&consensus);
        assert_eq!(result.sample, "sample1");
        assert_eq!(result.target, "target1");
        assert_eq!(result.variant_type, VariantType::Substitution);
        // rVAF should be updated to consensus value
        assert!((result.rvaf - 0.12).abs() < 1e-6);
        // info field should contain consensus metadata
        assert!(result.info.contains("multi_k="));
        assert!(result.info.contains("tier=Tier1"));
        assert!(result.info.contains("score=1.000"));
    }

    // -- types_compatible --

    #[test]
    fn test_types_compatible_same() {
        assert!(types_compatible(VariantType::Substitution, VariantType::Substitution));
        assert!(types_compatible(VariantType::Insertion, VariantType::Insertion));
    }

    #[test]
    fn test_types_compatible_indel_family() {
        assert!(types_compatible(VariantType::Insertion, VariantType::Deletion));
        assert!(types_compatible(VariantType::Insertion, VariantType::Indel));
        assert!(types_compatible(VariantType::Deletion, VariantType::Indel));
    }

    #[test]
    fn test_types_not_compatible_across_families() {
        assert!(!types_compatible(VariantType::Substitution, VariantType::Insertion));
        assert!(!types_compatible(VariantType::Substitution, VariantType::Deletion));
        assert!(!types_compatible(VariantType::Itd, VariantType::Substitution));
    }

    // -- Edge case: single k value --

    #[test]
    fn test_single_k_value_tier1() {
        let call = make_call(
            "s", "t", VariantType::Substitution,
            "5:A/T:5", 0.10, 100.0, 50,
            "A", "T", "ACGTACGT",
        );

        let calls_per_k = vec![(31, vec![call])];
        let result = merge_multi_k(&calls_per_k, 1);
        assert_eq!(result.len(), 1);
        // With total_k_values=1 and detected=1, this should be Tier1.
        assert_eq!(result[0].tier, ConfidenceTier::Tier1);
    }

    // -- Multiple variants from same k value --

    #[test]
    fn test_multiple_variants_same_k_stay_separate() {
        let snv1 = make_call(
            "s", "t", VariantType::Substitution,
            "3:A/T:3", 0.10, 100.0, 50,
            "A", "T", "ACGTACGT",
        );
        let snv2 = make_call(
            "s", "t", VariantType::Substitution,
            "5:C/G:5", 0.08, 80.0, 40,
            "C", "G", "ACGTACGT",
        );

        let calls_per_k = vec![(31, vec![snv1, snv2])];
        let result = merge_multi_k(&calls_per_k, 1);
        assert_eq!(result.len(), 2, "different variants should remain separate");
    }

    // -- Consensus score calculation --

    #[test]
    fn test_consensus_score_partial() {
        // For a Deletion at k=21 and k=31 (out of k=21,31,41):
        // weights: k21=1.0, k31=0.9, k41=0.7
        // actual = 1.0 + 0.9 = 1.9
        // max = 1.0 + 0.9 + 0.7 = 2.6
        // score = 1.9 / 2.6 ≈ 0.7308
        let call_k21 = make_call(
            "s", "t", VariantType::Deletion,
            "5:AT/A:6", 0.15, 90.0, 45,
            "AT", "A", "ACGTATCGT",
        );
        let call_k31 = make_call(
            "s", "t", VariantType::Deletion,
            "5:AT/A:6", 0.18, 110.0, 55,
            "AT", "A", "ACGTATCGT",
        );

        let calls_per_k = vec![
            (21, vec![call_k21]),
            (31, vec![call_k31]),
            (41, vec![]),
        ];

        let result = merge_multi_k(&calls_per_k, 3);
        assert_eq!(result.len(), 1);

        let expected_score = 1.9 / 2.6;
        assert!(
            (result[0].consensus_score - expected_score).abs() < 1e-4,
            "expected ~{:.4}, got {:.4}",
            expected_score,
            result[0].consensus_score,
        );
    }
}
