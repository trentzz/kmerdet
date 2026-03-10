# INDEL Filtering: The Representation Problem

## 1. The Core Problem

INDEL variant representation is inherently ambiguous. The same biological event -- a deletion or insertion of nucleotides -- can be described by multiple equivalent coordinate/allele combinations in VCF format. This creates a critical matching problem for tumor-informed filtering: a variant detected by k-mer walking may be represented differently than the same variant in the reference file, causing the filter to miss a true match.

### Concrete Example: Deletion in a Repeat

Consider a deletion of one "A" from a homopolymer run:

```
Reference:  5'-...GCTAAAACTT...-3'
            Position: 100 101 102 103 104 105 106 107 108 109
Mutant:     5'-...GCT-AAACTT...-3'  (one A deleted)
```

This deletion can be represented four equivalent ways in VCF:

| Representation | POS | REF | ALT | Description |
|----------------|-----|-----|-----|-------------|
| 1 | 100 | TA | T | Delete A at position 101 |
| 2 | 101 | AA | A | Delete A at position 102 |
| 3 | 102 | AA | A | Delete A at position 103 |
| 4 | 103 | AA | A | Delete A at position 104 |

All four describe the identical biological event. VCF convention requires position to be 1-based and to include the preceding reference base. The problem is that different tools may choose different representations:

- **k-mer walking** reports whichever representation the graph pathfinding happened to produce, which depends on the walking direction and which k-mer first diverged from the reference
- **The reference file** may use whichever representation the original variant caller produced (GATK, VarScan, manual annotation, etc.)

If km finds representation 3 and the reference file uses representation 1, coordinate-based filtering (reference mode) reports "Not Found" -- a false negative in the filter.

### Scale of the Problem

From the thesis validation study:
- INDEL sensitivity was significantly lower than SNV sensitivity (77% for SNVs)
- A portion of INDEL false negatives were attributable to representation mismatches, not detection failures
- The `use-alt` mode in kmtools was developed specifically to address this problem, but it was only available in the development branch and was not used in the initial validation

SNVs do not have this problem: a single-base substitution has exactly one representation (one POS, one REF base, one ALT base). INDELs are the only variant class with ambiguous coordinate representation.

## 2. Left-Alignment Normalization

### The Standard Solution

Left-alignment normalization shifts INDEL representations to the leftmost possible position while maintaining allele length. Combined with parsimony (using the shortest possible REF and ALT alleles), this produces a unique canonical representation for every INDEL.

A variant is considered **normalized** if and only if:
1. **It is parsimonious:** The REF and ALT alleles share no common suffix, and the shortest possible alleles are used (Tan et al., 2015).
2. **It is left-aligned:** It is not possible to shift the variant position one base to the left while keeping allele lengths constant.

### Algorithm (Tan et al., 2015)

```
function normalize(chrom, pos, ref, alt):
    # Step 1: Trim common suffix
    while ref[-1] == alt[-1] and len(ref) > 1 and len(alt) > 1:
        ref = ref[:-1]
        alt = alt[:-1]

    # Step 2: Trim common prefix (keep at least 1 base for anchor)
    while ref[0] == alt[0] and len(ref) > 1 and len(alt) > 1:
        ref = ref[1:]
        alt = alt[1:]
        pos += 1

    # Step 3: Left-align
    while pos > 1:
        preceding_base = reference_genome[chrom][pos - 2]  # 0-based
        if ref[-1] == preceding_base:
            # Can shift left
            ref = preceding_base + ref[:-1]
            alt = preceding_base + alt[:-1]
            pos -= 1
        else:
            break

    # Step 4: Re-trim (left-alignment may introduce new common prefix/suffix)
    # ... repeat steps 1-2

    return (chrom, pos, ref, alt)
```

### Available Tools

- **bcftools norm** (`bcftools norm -f ref.fa -m -both`): Left-aligns and splits multiallelic sites. The most widely used normalization tool.
- **vt normalize** (`vt normalize -r ref.fa`): Fast normalizer from the vt toolkit, specifically designed for INDEL normalization.
- **GATK LeftAlignAndTrimVariants**: Java-based normalizer integrated into the GATK suite.

All three produce identical normalized representations for well-formed VCF records.

### Implementation in kmerdet

kmerdet should normalize INDELs at two points:

1. **During detection (output normalization):** After classifying a variant as an insertion or deletion, normalize its coordinates before writing to output.
2. **During filtering (input normalization):** Normalize both the detected variant and the expected variant from the reference file before comparing.

```rust
pub fn normalize_indel(
    chrom: &str,
    pos: u64,
    ref_allele: &str,
    alt_allele: &str,
    reference: &ReferenceGenome,
) -> NormalizedVariant {
    let mut pos = pos;
    let mut ref_bytes = ref_allele.as_bytes().to_vec();
    let mut alt_bytes = alt_allele.as_bytes().to_vec();

    // Step 1: Trim common suffix
    while ref_bytes.len() > 1 && alt_bytes.len() > 1
        && ref_bytes.last() == alt_bytes.last()
    {
        ref_bytes.pop();
        alt_bytes.pop();
    }

    // Step 2: Left-align
    loop {
        if pos <= 1 { break; }
        let preceding = reference.base_at(chrom, pos - 2); // 0-based
        if ref_bytes.last() == Some(&preceding) {
            ref_bytes.pop();
            ref_bytes.insert(0, preceding);
            alt_bytes.pop();
            alt_bytes.insert(0, preceding);
            pos -= 1;
        } else {
            break;
        }
    }

    // Step 3: Trim common prefix (keep anchor)
    while ref_bytes.len() > 1 && alt_bytes.len() > 1
        && ref_bytes[0] == alt_bytes[0]
    {
        ref_bytes.remove(0);
        alt_bytes.remove(0);
        pos += 1;
    }

    NormalizedVariant {
        chrom: chrom.to_string(),
        pos,
        ref_allele: String::from_utf8(ref_bytes).unwrap(),
        alt_allele: String::from_utf8(alt_bytes).unwrap(),
    }
}
```

**Reference genome requirement:** Left-alignment requires access to the reference genome sequence to know what base precedes the variant. For kmerdet, this can come from:
- The target FASTA file (if the variant is within the target sequence)
- A full reference genome file (for variants near target boundaries)

For targets with sufficient flanking sequence (the standard `--flank 35` provides 35 bp on each side), the target FASTA itself is sufficient. This avoids requiring a full reference genome as a dependency.

## 3. ALT Sequence Comparison

### The Generalized Solution

Instead of matching by coordinates (CHROM + POS + REF + ALT), compare the full alternative sequence produced by k-mer walking against the expected variant's alternative sequence. This is what kmtools' `use-alt` mode implements.

**How it works:**

1. k-mer walking produces a full variant path: a sequence of k-mers that, when overlapped, spell out the variant sequence
2. The reference file specifies the expected variant, from which the expected alternative sequence can be reconstructed
3. Matching: does the detected alternative sequence contain (or equal) the expected alternative sequence?

```rust
pub fn alt_sequence_match(
    detected_alt_seq: &str,
    expected_alt_seq: &str,
    min_overlap: usize,  // minimum shared bases for a match
) -> bool {
    // Exact match
    if detected_alt_seq == expected_alt_seq {
        return true;
    }
    // Substring match (detected contains expected or vice versa)
    if detected_alt_seq.contains(expected_alt_seq)
        || expected_alt_seq.contains(detected_alt_seq)
    {
        return true;
    }
    // Overlap match (for partial detections)
    let overlap = longest_common_substring(detected_alt_seq, expected_alt_seq);
    overlap >= min_overlap
}
```

### Advantages Over Coordinate Matching

| Aspect | Coordinate Matching | ALT Sequence Matching |
|--------|-------------------|-----------------------|
| INDEL representation | Sensitive to representation | Representation-agnostic |
| Complex INDELs | Fails for delins with different decomposition | Matches on full sequence |
| MNVs | Requires exact position match | Matches on altered bases |
| Computational cost | O(1) hash lookup | O(n*m) string comparison |
| False match risk | Very low | Moderate (subsequence may match by chance) |

### When ALT Sequence Matching Fails

ALT sequence matching is not a panacea:

- **Short INDELs near repeats:** A 1-bp insertion in a long repeat creates an alternative sequence almost identical to the reference. The inserted base is indistinguishable from the surrounding context, making substring matching unreliable.
- **Multiple nearby variants:** If two variants are close together, the alternative sequence from one variant may partially overlap with the other, causing incorrect matches.
- **Partial detection:** If k-mer walking detects only part of a complex variant (e.g., one end of a delins), the alternative sequence is incomplete and may not match the expected full sequence.

## 4. Normalization as Default Behavior

### The Case for Always Normalizing

The current design in `src/filter/reference_mode.rs` performs coordinate matching without normalization. If the detected variant and expected variant have different representations, the match fails silently. The user receives "Not Found" with no indication that the variant was actually detected but in a different representation.

**Proposed change:** Always normalize INDEL coordinates before comparison in reference mode:

```rust
pub fn find_match_normalized<'a>(
    expected: &ExpectedVariant,
    calls: &'a [VariantCall],
    reference: &ReferenceGenome,
) -> Option<&'a VariantCall> {
    let norm_expected = normalize_indel(
        &expected.chrom, expected.pos,
        &expected.ref_allele, &expected.alt_allele,
        reference,
    );

    calls.iter().find(|call| {
        let norm_call = normalize_indel(
            &call.chrom, call.pos,
            &call.ref_allele, &call.alt_allele,
            reference,
        );
        norm_expected.chrom == norm_call.chrom
            && norm_expected.pos == norm_call.pos
            && norm_expected.ref_allele == norm_call.ref_allele
            && norm_expected.alt_allele == norm_call.alt_allele
    })
}
```

### Normalization + ALT Fallback Strategy

The recommended matching strategy uses normalization first, then falls back to ALT sequence comparison if normalization does not produce a match:

```
Step 1: Normalize both detected and expected variants
Step 2: Attempt coordinate match on normalized representations
Step 3: If no match, attempt ALT sequence comparison
Step 4: Report match source (coordinate vs. alt-sequence) in output
```

This maximizes sensitivity while preserving the precision of coordinate matching. The fallback to ALT matching handles cases where normalization alone is insufficient (e.g., complex INDELs where the k-mer-derived representation differs structurally from the expected representation).

### Backward Compatibility

For backward compatibility with the thesis pipeline, normalization should be enabled by default but disableable:

```
kmerdet filter --no-normalize  # disable normalization (legacy behavior)
kmerdet filter                 # normalization enabled (new default)
kmerdet filter --use-alt       # ALT sequence matching (overrides coordinate matching)
kmerdet filter --normalize-then-alt  # normalization first, ALT fallback (recommended)
```

## 5. Handling Equivalent INDEL Representations in K-mer Space

### Multiple Paths for the Same INDEL

K-mer walking may discover the same INDEL through different graph paths, producing multiple variant calls that represent the same biological event. For example, a deletion in a repeat region may be found by branching off the reference at different positions within the repeat.

**Example:** Deletion of one "T" from "TTTTT" at positions 200-204:

```
Path 1: branch at position 200 → skip T → rejoin at 201
Path 2: branch at position 201 → skip T → rejoin at 202
Path 3: branch at position 202 → skip T → rejoin at 203
```

All three paths represent the same deletion. Without deduplication, the variant is reported three times with slightly different coordinates and potentially different rVAF estimates.

### Deduplication by Normalization

Normalize all detected INDELs before reporting:

```rust
pub fn deduplicate_indels(calls: &mut Vec<VariantCall>, reference: &ReferenceGenome) {
    // Normalize all calls
    for call in calls.iter_mut() {
        if call.variant_type == "Insertion" || call.variant_type == "Deletion" {
            let norm = normalize_indel(
                &call.chrom, call.pos,
                &call.ref_allele, &call.alt_allele,
                reference,
            );
            call.pos = norm.pos;
            call.ref_allele = norm.ref_allele;
            call.alt_allele = norm.alt_allele;
        }
    }

    // Deduplicate by (chrom, pos, ref, alt, type)
    calls.sort_by(|a, b| {
        (&a.chrom, a.pos, &a.ref_allele, &a.alt_allele)
            .cmp(&(&b.chrom, b.pos, &b.ref_allele, &b.alt_allele))
    });
    calls.dedup_by(|a, b| {
        a.chrom == b.chrom && a.pos == b.pos
            && a.ref_allele == b.ref_allele
            && a.alt_allele == b.alt_allele
            && a.variant_type == b.variant_type
    });
}
```

When deduplicating, keep the call with the highest rVAF (or expression), as it likely represents the best-supported observation.

### Canonical INDEL Representation

Define a canonical form for every INDEL:

```
Canonical = (chromosome, left_aligned_position, parsimonious_ref, parsimonious_alt)
```

This is the output of the normalization algorithm. All internal operations (deduplication, filtering, merging, PoN lookup) should use the canonical form. The original (un-normalized) representation can be retained in a separate field for debugging.

## 6. Complex INDEL Challenges

### Delins (Deletion-Insertion)

A delins is a simultaneous deletion and insertion at the same position. For example:

```
Reference:  ACGTACGT
Mutant:     AC--TTCGT  (delete GT, insert TT → net: G>T substitution spanning 2 bases)
```

This could be represented as:
- **MNV:** pos=3, REF=GT, ALT=TT (two substitutions)
- **Delins:** pos=3, REF=GTA, ALT=TT (delete 3 bases, insert 2)
- **Complex:** pos=3, REF=GTAC, ALT=TTAC (include trailing context)

K-mer walking does not know which interpretation is correct. The `diff_path_without_overlap` classification logic attempts to determine the variant type from the path comparison, but for delins the distinction between an MNV and a delins is inherently ambiguous in k-mer space.

**Recommendation:** Report delins as `Complex` type with both the deleted and inserted sequences. Let downstream annotation tools (VEP, ANNOVAR) determine the biological interpretation.

### MNVs Near INDELs

When a substitution and an INDEL are within k bases of each other, they may appear as a single complex variant in the k-mer graph. The graph has a single branching-merging path that encompasses both events.

**Phase matters:** If both variants are on the same allele (cis), they produce a single alternative path. If they are on different alleles (trans), they produce two separate alternative paths. The k-mer graph topology distinguishes these cases:

```
Cis (same allele):
  Ref: ─────────────────
       └── MNV+INDEL ──┘  (single alt path)

Trans (different alleles):
  Ref: ─────────────────
       ├── MNV only ───┤  (alt path 1)
       └── INDEL only ─┘  (alt path 2)
```

If the cluster mode correctly identifies both paths, the NNLS quantification separates the contributions. However, at very low VAF, one or both paths may not have sufficient support for reliable decomposition.

### Tandem Duplications

Tandem duplications (TDs) are a special case where the inserted sequence is a copy of the adjacent reference sequence. In VCF:

```
POS=100, REF=A, ALT=AACGT  (insert ACGT, which is a copy of positions 101-104)
```

This can be represented either as:
- **Insertion at position 100:** POS=100, REF=A, ALT=AACGT
- **Duplication call:** DUP at 101-104

K-mer walking sees a tandem duplication as an insertion, but the classification logic in km can detect it as an ITD (Internal Tandem Duplication) when `start == end_ref_overlap` -- meaning the variant path retraces a segment of the reference.

**Normalization consideration:** Left-alignment of a tandem duplication shifts the insertion to the leftmost position in the tandem repeat unit. If the duplicated segment is repeated multiple times, left-alignment places the insertion at the first copy. This is consistent with VCF convention but may differ from clinical reporting conventions (e.g., FLT3-ITD is typically reported at the exon 14-15 junction regardless of exact insertion position).

### Annotation-Dependent Interpretation

Different annotation frameworks interpret the same variant differently:

| Framework | Interpretation of delins at pos 100, REF=ACG, ALT=TT |
|-----------|------------------------------------------------------|
| VCF literal | Complex INDEL: delete ACG, insert TT |
| HGVS | c.100_102delinsTT |
| Protein effect | Depends on reading frame and amino acid impact |
| Clinical report | May simplify to "deletion at codon X" |

kmerdet should report variants in a framework-neutral format (VCF-style coordinates) and leave interpretation to downstream annotation. The output should include the raw reference and alternative sequences to enable any interpretation framework.

## 7. Testing INDEL Normalization

### Test Cases

A comprehensive test suite for INDEL normalization should include:

```rust
#[cfg(test)]
mod tests {
    // Simple deletion in homopolymer
    // Input: pos=103, ref=AA, alt=A in GCTAAAACTT
    // Expected: pos=101, ref=TA, alt=T (left-aligned)

    // Insertion in homopolymer
    // Input: pos=105, ref=A, alt=AA in GCTAAAACTT
    // Expected: pos=100, ref=G, alt=GA (left-aligned)

    // Complex delins (no shift possible)
    // Input: pos=100, ref=ACG, alt=TT
    // Expected: pos=100, ref=ACG, alt=TT (already normalized)

    // Deletion at start of sequence (cannot left-align further)
    // Input: pos=1, ref=AT, alt=A
    // Expected: pos=1, ref=AT, alt=A (already at leftmost position)

    // Multi-base deletion in tandem repeat
    // Input: pos=102, ref=ATAT, alt=AT in ...ATATATATATAT...
    // Expected: pos=1, ref=ATAT, alt=AT (shifted to start of repeat)

    // SNV (should not be modified)
    // Input: pos=100, ref=A, alt=T
    // Expected: pos=100, ref=A, alt=T (unchanged)
}
```

### Validation Against bcftools

For integration testing, normalize a set of VCF records with both kmerdet and `bcftools norm`, and verify identical output:

```bash
# Create test VCF with known INDELs
# Run bcftools norm
bcftools norm -f reference.fa test_indels.vcf -o bcftools_normalized.vcf

# Run kmerdet normalization
kmerdet normalize test_indels.vcf -r reference.fa -o kmerdet_normalized.vcf

# Compare (ignoring header differences)
diff <(grep -v '^#' bcftools_normalized.vcf) <(grep -v '^#' kmerdet_normalized.vcf)
```

Any differences indicate a bug in kmerdet's normalization implementation.

### Regression Tests from Thesis Data

The thesis 10-patient validation identified specific INDEL matching failures. These should be preserved as regression tests:

- Patient-specific INDELs where reference mode failed but the variant was genuinely detected
- Known INDEL representations from the reference file and the corresponding km output representation
- The expected normalized form for each

These tests ensure that the normalization fix actually resolves the matching failures observed in the thesis.

## References

- Tan et al. (2015) -- Unified representation of genetic variants (Bioinformatics; defines parsimony and left-alignment)
- Danecek et al. (2021) -- Twelve years of SAMtools and BCFtools (includes bcftools norm specification)
- Li (2014) -- vt: a tool set for short variant discovery (vt normalize algorithm)
- Van der Auwera & O'Connor (2020) -- Genomics in the Cloud: GATK (LeftAlignAndTrimVariants)
- VCF specification v4.3 -- Variant Call Format (representation conventions)
- Audoux et al. (2017) -- km: diff_path_without_overlap classification logic
- UPS-indel (Hasan et al., 2017) -- Universal Positioning System for INDELs (alternative normalization approach)
