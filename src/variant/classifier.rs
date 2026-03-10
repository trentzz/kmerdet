/// Variant classification from reference/alternative path comparison.
///
/// Implements `diff_path_without_overlap` from km:
/// 1. Left-to-right scan for first divergence point (i)
/// 2. Right-to-left scan for last divergence point (j_ref, j_var)
/// 3. Classify based on relative positions and lengths

use super::VariantType;
use crate::sequence::path::KmerPath;

/// Result of comparing a reference path to an alternative path.
#[derive(Debug)]
pub struct Classification {
    pub variant_type: VariantType,
    /// Human-readable variant name (e.g., "41:A/T:41" for SNV).
    pub variant_name: String,
    /// Reference allele at the variant position.
    pub ref_allele: String,
    /// Alternative allele at the variant position.
    pub alt_allele: String,
    /// Start position of the variant within the target.
    pub start: usize,
    /// End position of the variant within the target.
    pub end: usize,
}

/// Classify the variant by comparing reference and alternative paths.
///
/// Implements the `diff_path_without_overlap` algorithm:
/// 1. Convert paths to DNA sequences
/// 2. Find first and last divergence points
/// 3. Extract ref/alt alleles
/// 4. Classify based on allele lengths and content
pub fn classify(ref_path: &KmerPath, alt_path: &KmerPath, _k: u8) -> Classification {
    let ref_seq: Vec<u8> = ref_path.to_sequence().into_bytes();
    let alt_seq: Vec<u8> = alt_path.to_sequence().into_bytes();

    // Step 2: Identical sequences → Reference
    if ref_seq == alt_seq {
        return Classification {
            variant_type: VariantType::Reference,
            variant_name: String::new(),
            ref_allele: String::new(),
            alt_allele: String::new(),
            start: 0,
            end: 0,
        };
    }

    // Step 3: Find first divergence point (left to right)
    let start = ref_seq
        .iter()
        .zip(alt_seq.iter())
        .position(|(r, a)| r != a)
        .unwrap_or(std::cmp::min(ref_seq.len(), alt_seq.len()));

    // Step 4: Find last divergence point (right to left)
    let mut j_ref = ref_seq.len().saturating_sub(1);
    let mut j_alt = alt_seq.len().saturating_sub(1);

    while j_ref > start && j_alt > start && ref_seq[j_ref] == alt_seq[j_alt] {
        j_ref -= 1;
        j_alt -= 1;
    }

    let end_ref = j_ref;
    let end_alt = j_alt;

    // Step 5: Extract alleles
    let ref_allele = if end_ref >= start && start < ref_seq.len() {
        String::from_utf8_lossy(&ref_seq[start..=end_ref.min(ref_seq.len() - 1)]).to_string()
    } else {
        String::new()
    };

    let alt_allele = if end_alt >= start && start < alt_seq.len() {
        String::from_utf8_lossy(&alt_seq[start..=end_alt.min(alt_seq.len() - 1)]).to_string()
    } else {
        String::new()
    };

    // Step 6: Classify
    let ref_len = ref_allele.len();
    let alt_len = alt_allele.len();

    let variant_type = if ref_len == alt_len {
        // Equal lengths: substitution (SNV or MNV)
        VariantType::Substitution
    } else if alt_len > ref_len {
        // Alt is longer: could be insertion, ITD, or indel
        if is_itd(&ref_allele, &alt_allele, &ref_seq, start) {
            VariantType::Itd
        } else if ref_len == 0 {
            // Pure insertion (no ref bases replaced)
            VariantType::Insertion
        } else {
            // Both have content but different lengths → complex indel or insertion
            // If ref_len == 1, this is effectively a pure insertion with an anchor base
            // Check if it's a simple insertion (ref is a prefix/suffix of alt)
            if alt_allele.starts_with(&ref_allele) || alt_allele.ends_with(&ref_allele) {
                VariantType::Insertion
            } else {
                VariantType::Indel
            }
        }
    } else {
        // Ref is longer: could be deletion or indel
        if alt_len == 0 {
            // Pure deletion
            VariantType::Deletion
        } else if ref_allele.starts_with(&alt_allele) || ref_allele.ends_with(&alt_allele) {
            VariantType::Deletion
        } else {
            VariantType::Indel
        }
    };

    // Step 7: Build variant name
    let end = std::cmp::max(end_ref, end_alt);
    let variant_name = format!("{}:{}/{}:{}", start, ref_allele, alt_allele, end);

    Classification {
        variant_type,
        variant_name,
        ref_allele,
        alt_allele,
        start,
        end,
    }
}

/// Check if the variant is an internal tandem duplication (ITD).
///
/// ITD occurs when the alt allele contains a tandem repeat of a segment from
/// the reference. We check if the extra bases in alt_allele are an exact copy
/// of adjacent reference sequence.
fn is_itd(ref_allele: &str, alt_allele: &str, ref_seq: &[u8], start: usize) -> bool {
    if alt_allele.len() <= ref_allele.len() {
        return false;
    }

    let alt_bytes = alt_allele.as_bytes();
    let extra_len = alt_bytes.len() - ref_allele.len();

    // ITD requires at least 2 duplicated bases; single-base duplications
    // are classified as simple insertions rather than tandem duplications.
    if extra_len < 2 {
        return false;
    }

    // Pure insertion (empty ref allele): check if inserted bases duplicate
    // adjacent reference sequence.
    if ref_allele.is_empty() {
        let ins_len = alt_bytes.len();

        // Check if inserted bases match ref sequence just before the insertion
        if start >= ins_len {
            let ref_segment = &ref_seq[start - ins_len..start];
            if alt_bytes == ref_segment {
                return true;
            }
        }
        // Check if inserted bases match ref sequence starting at insertion point
        if start + ins_len <= ref_seq.len() {
            let ref_segment = &ref_seq[start..start + ins_len];
            if alt_bytes == ref_segment {
                return true;
            }
        }
        return false;
    }

    // Non-empty ref allele: the alt allele is longer than the ref allele.
    // Check if the entire alt allele can be explained as the ref allele
    // plus a tandem duplicate of adjacent reference sequence.

    // Strategy: the alt allele occupies the variant region. Extract the
    // "extra" bases (the ones beyond the ref allele length). Check whether
    // these extra bases match a segment of the reference sequence adjacent
    // to the variant region, which would indicate tandem duplication.

    let ref_bytes = ref_allele.as_bytes();

    // Check if alt = ref + extra, where extra duplicates reference at `start`
    if alt_allele.starts_with(ref_allele) {
        let extra = &alt_bytes[ref_bytes.len()..];
        if start + extra_len <= ref_seq.len() && extra == &ref_seq[start..start + extra_len] {
            return true;
        }
    }

    // Check if alt = extra + ref, where extra duplicates reference before variant end
    if alt_allele.ends_with(ref_allele) {
        let extra = &alt_bytes[..extra_len];
        // Extra matches ref just before start
        if start >= extra_len && extra == &ref_seq[start - extra_len..start] {
            return true;
        }
        // Extra matches ref just after the variant region
        let ref_end = start + ref_bytes.len();
        if ref_end + extra_len <= ref_seq.len() && extra == &ref_seq[ref_end..ref_end + extra_len]
        {
            return true;
        }
    }

    // General check: the alt allele might contain the ref allele in the middle
    // with extra bases that duplicate adjacent reference. Try all split points
    // where we can decompose alt_allele into a portion matching ref and an
    // "extra" portion matching adjacent reference.
    //
    // The key insight: the entire alt_allele region corresponds to ref_seq[start..=end_ref].
    // If the alt_allele is a tandem duplication, then the alt_allele should consist of
    // reference bases from a wider window around the variant position.
    //
    // Check if alt_allele matches a contiguous segment of the reference that
    // includes the variant region plus adjacent bases (i.e., the alt is simply
    // a copy of a longer reference window).
    // Check if alt matches ref_seq[start - offset .. start - offset + alt_len] for some offset
    for offset in 0..=extra_len {
        if start >= offset {
            let window_start = start - offset;
            let window_end = window_start + alt_bytes.len();
            if window_end <= ref_seq.len() && alt_bytes == &ref_seq[window_start..window_end] {
                return true;
            }
        }
    }

    false
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: build a KmerPath that produces the given sequence for the given k.
    fn make_path(seq: &str, k: usize, is_reference: bool) -> KmerPath {
        assert!(seq.len() >= k, "sequence must be at least k bases long");
        let kmers: Vec<String> = (0..=seq.len() - k)
            .map(|i| seq[i..i + k].to_string())
            .collect();
        KmerPath {
            kmers,
            is_reference,
        }
    }

    #[test]
    fn test_classify_reference() {
        let ref_path = make_path("ACGTACGT", 4, true);
        let alt_path = make_path("ACGTACGT", 4, false);
        let result = classify(&ref_path, &alt_path, 4);
        assert_eq!(result.variant_type, VariantType::Reference);
    }

    #[test]
    fn test_classify_snv() {
        // A→T at position 4
        let ref_path = make_path("ACGTACGT", 4, true);
        let alt_path = make_path("ACGTTCGT", 4, false);
        let result = classify(&ref_path, &alt_path, 4);
        assert_eq!(result.variant_type, VariantType::Substitution);
        assert_eq!(result.ref_allele, "A");
        assert_eq!(result.alt_allele, "T");
        assert_eq!(result.start, 4);
    }

    #[test]
    fn test_classify_mnv() {
        // ACG→TCA at positions 2-4
        let ref_path = make_path("TTACGTTT", 4, true);
        let alt_path = make_path("TTTCATTT", 4, false);
        let result = classify(&ref_path, &alt_path, 4);
        assert_eq!(result.variant_type, VariantType::Substitution);
        assert_eq!(result.ref_allele, "ACG");
        assert_eq!(result.alt_allele, "TCA");
        assert_eq!(result.start, 2);
    }

    #[test]
    fn test_classify_insertion_single_base() {
        // A inserted between pos 4 and 5
        let ref_path = make_path("ACGTACGT", 4, true);
        let alt_path = make_path("ACGTAACGT", 4, false);
        let result = classify(&ref_path, &alt_path, 4);
        assert!(
            result.variant_type == VariantType::Insertion,
            "expected Insertion, got {:?}",
            result.variant_type
        );
    }

    #[test]
    fn test_classify_insertion_multi_base() {
        // GGG inserted
        let ref_path = make_path("ACGTACGT", 4, true);
        let alt_path = make_path("ACGTGGGACGT", 4, false);
        let result = classify(&ref_path, &alt_path, 4);
        assert!(
            result.variant_type == VariantType::Insertion,
            "expected Insertion, got {:?}",
            result.variant_type
        );
        assert!(result.alt_allele.len() > result.ref_allele.len());
    }

    #[test]
    fn test_classify_deletion_single_base() {
        // T deleted at pos 3
        let ref_path = make_path("ACGTACGT", 4, true);
        let alt_path = make_path("ACGACGT", 4, false);
        let result = classify(&ref_path, &alt_path, 4);
        assert!(
            result.variant_type == VariantType::Deletion,
            "expected Deletion, got {:?}",
            result.variant_type
        );
    }

    #[test]
    fn test_classify_deletion_multi_base() {
        // GT deleted
        let ref_path = make_path("ACGTACGT", 4, true);
        let alt_path = make_path("ACACGT", 4, false);
        let result = classify(&ref_path, &alt_path, 4);
        assert!(
            result.variant_type == VariantType::Deletion,
            "expected Deletion, got {:?}",
            result.variant_type
        );
        assert!(result.ref_allele.len() > result.alt_allele.len());
    }

    #[test]
    fn test_classify_complex_indel() {
        // Complex indel: different bases AND different lengths
        // ref: ACGTACGT (8)  → alt: ACTTACGT (8): same len, would be substitution
        // Need different lengths. Use:
        // ref: ACGTAACGT (9) → alt: ACTCCGT (7): GTA→TC (3 ref bases → 2 different alt bases)
        let ref_path = make_path("ACGTAACGT", 4, true);
        let alt_path = make_path("ACTCCGT", 4, false);
        let result = classify(&ref_path, &alt_path, 4);
        assert_eq!(
            result.variant_type,
            VariantType::Indel,
            "expected Indel, got {:?} (ref='{}', alt='{}')",
            result.variant_type,
            result.ref_allele,
            result.alt_allele,
        );
    }

    #[test]
    fn test_classify_itd() {
        // ACGT duplicated: ref ACGTACGT → alt ACGTACGTACGT
        let ref_path = make_path("ACGTACGT", 4, true);
        let alt_path = make_path("ACGTACGTACGT", 4, false);
        let result = classify(&ref_path, &alt_path, 4);
        assert_eq!(
            result.variant_type,
            VariantType::Itd,
            "expected ITD, got {:?}",
            result.variant_type
        );
    }

    #[test]
    fn test_classify_variant_name_format() {
        // SNV: should produce "start:ref/alt:end" format
        let ref_path = make_path("ACGTACGT", 4, true);
        let alt_path = make_path("ACGTTCGT", 4, false);
        let result = classify(&ref_path, &alt_path, 4);

        // Verify format: "start:ref_allele/alt_allele:end"
        let parts: Vec<&str> = result.variant_name.splitn(2, ':').collect();
        assert!(parts.len() >= 2, "variant_name should contain ':'");

        // Parse: "start:ref/alt:end"
        let name = &result.variant_name;
        assert!(
            name.contains('/'),
            "variant_name should contain '/' separator"
        );
        assert!(
            name.starts_with(&result.start.to_string()),
            "variant_name should start with start position"
        );
        assert!(
            name.ends_with(&result.end.to_string()),
            "variant_name should end with end position"
        );
    }

    #[test]
    fn test_classify_alleles_correct() {
        // Verify that ref_allele and alt_allele are exactly the changed bases
        // SNV at position 4: A→T
        let ref_path = make_path("ACGTACGT", 4, true);
        let alt_path = make_path("ACGTTCGT", 4, false);
        let result = classify(&ref_path, &alt_path, 4);
        assert_eq!(result.ref_allele, "A", "ref_allele should be just 'A'");
        assert_eq!(result.alt_allele, "T", "alt_allele should be just 'T'");
        assert_eq!(result.start, 4);
        assert_eq!(result.end, 4);

        // Deletion: ACGTACGT → ACGACGT (T at pos 3 deleted)
        let ref_path2 = make_path("ACGTACGT", 4, true);
        let alt_path2 = make_path("ACGACGT", 4, false);
        let result2 = classify(&ref_path2, &alt_path2, 4);
        // The ref allele should contain the deleted base(s)
        assert!(
            !result2.ref_allele.is_empty(),
            "ref_allele should not be empty for deletion"
        );
    }
}
