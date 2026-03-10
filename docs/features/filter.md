# Feature: `filter` -- Tumor-Informed Filtering

## What It Does

The `filter` subcommand applies tumor-informed filtering to detection results, retaining
only variants that match a patient's known mutation profile. It operates in two modes:
reference mode (coordinate-based matching) and use-alt mode (sequence-based matching).
INDEL normalization is applied by default before matching to handle equivalent variant
representations.

```
kmerdet filter -i results.tsv -r reference-info.tsv --mode reference -o filtered.tsv
kmerdet filter -i results.tsv -r reference-info.tsv --mode use-alt -o filtered.tsv
```

### Input Specification

| Input | Flag | Format | Required |
|-------|------|--------|----------|
| Detection results | `-i` / `--input` | TSV/CSV from `detect` or stdin | Yes |
| Reference info | `-r` / `--reference` | TSV with CHROM/POS/REF/ALT/TYPE columns | Yes (reference mode) |
| Additional references | `-r` (repeated) | Multiple reference files | No |
| Panel of normals | `--pon` | Binary PoN database | No |

### Output Specification

Same format as input, with an additional `FILTER` column indicating PASS or the reason
for filtering. The `--keep-filtered` flag retains all rows with filter annotations;
by default, only PASS rows are written.

---

## Why It Matters

Tumor-informed filtering is the final safeguard between raw variant calls and clinical
reporting. The thesis validation achieved zero false positives specifically because of
the tumor-informed filtering step: only variants matching the patient's known somatic
mutation profile were reported. Without filtering, the detection step would produce
variant calls for any non-reference path in the k-mer graph, including germline variants,
sequencing artifacts, and CHIP mutations.

For MRD monitoring, the clinical question is specific: "Is this patient's known tumor
mutation detectable in the current blood sample?" The filter answers this by comparing
detected variants against the patient's tumor mutation list. A detected variant that
does not match any known tumor mutation is irrelevant to the MRD question and is
filtered out.

The INDEL representation problem, documented extensively in indel-filtering.md, is a
primary source of false negatives in filtering. The same biological INDEL can be
represented with different coordinates depending on the tool that called it. The thesis
found that a portion of apparent INDEL false negatives were actually detection successes
that failed at the filtering step due to representation mismatches. Normalization
addresses this directly.

---

## Reference Mode

Reference mode matches variants by genomic coordinates and allele identity. All five
fields must match for a positive hit:

| Field | Comparison | Notes |
|-------|-----------|-------|
| CHROM | Exact string match | Must handle chr-prefix variation (chr17 vs 17) |
| POS | Exact integer match | After normalization |
| REF | Exact string match | After normalization |
| ALT | Exact string match | After normalization |
| TYPE | Exact string match | Substitution, Insertion, Deletion, ITD, Complex |

### Reference Information File Format

```tsv
CHROM	POS	REF	ALT	TYPE	GENE	ANNOTATION
chr17	7577120	C	T	Substitution	TP53	R248W
chr5	170837543	TCTG	T	Deletion	NPM1	W288fs
chr13	28018505	A	AACGT	Insertion	FLT3	ITD_exon14
```

The GENE and ANNOTATION columns are optional metadata passed through to the output
for clinical context.

### Chromosome Name Normalization

To handle inconsistency between reference files and detection output:
- Strip leading "chr" prefix and compare numerically where possible
- `chr17` matches `17`, `chrX` matches `X`, `chrM` matches `MT`
- Configurable via `--chr-style` (ucsc | ensembl | auto)

---

## Use-Alt Mode

Use-alt mode matches variants by the full alternative sequence rather than genomic
coordinates. This is the robust fallback for cases where coordinate normalization
is insufficient, particularly complex INDELs where km and the reference caller
produce structurally different representations of the same event.

### Matching Logic

```
1. Extract the alternative sequence from the detected variant path
2. Extract the expected alternative sequence from the reference info
3. Match if:
   a. Exact string equality, OR
   b. One sequence contains the other as a substring, OR
   c. Longest common substring >= min_overlap (default: k-1 = 30 bases)
```

Use-alt mode is more permissive than reference mode and should be used when reference
mode produces unexpected false negatives for INDELs. It is selected via `--mode use-alt`.

---

## INDEL Normalization

INDEL normalization is applied by default before coordinate matching in reference mode.
The algorithm follows Tan et al. (2015):

1. **Trim common suffix** from REF and ALT alleles
2. **Left-align** by shifting the variant leftward as long as the trailing base
   matches the preceding reference base
3. **Trim common prefix** (keeping one anchor base per VCF convention)

Both the detected variant and the expected variant from the reference file are
normalized before comparison. This ensures that equivalent representations produce
the same canonical form.

### Why Normalization Is Critical

The indel-filtering.md research document details the core problem: a single A deletion
in a poly-A run of length 4 has four equivalent VCF representations, differing only in
POS. K-mer walking may produce any of these depending on which k-mer first diverged
from the reference. The reference file may use whichever representation the original
caller produced. Without normalization, matching fails silently.

### Reference Sequence for Left-Alignment

Left-alignment requires knowing the preceding reference base. For targets with standard
flanking (35 bp on each side), the target FASTA itself provides sufficient context. When
the variant is near the target boundary and more context is needed, a reference genome
file can be supplied via `--ref-genome`.

### Disabling Normalization

For backward compatibility with the thesis pipeline, normalization can be disabled:

```
kmerdet filter --no-normalize    # Legacy coordinate matching without normalization
```

---

## Additional Filters

Beyond tumor-informed matching, the filter subcommand supports threshold-based filters
that can be combined with either mode:

| Filter | Flag | Default | Description |
|--------|------|---------|-------------|
| Minimum rVAF | `--min-vaf` | 0.001 | Remove calls below this VAF |
| Minimum expression | `--min-expression` | 5.0 | Minimum path coefficient |
| Minimum coverage | `--min-coverage` | 2 | Minimum k-mer count along path |
| Exclude types | `--exclude-type` | ["Reference"] | Variant types to remove |
| Maximum strand bias | `--max-strand-bias` | None | Filter calls with extreme strand bias |
| Minimum confidence | `--min-confidence` | None | Phred-scaled quality threshold |

Filters are applied in sequence after tumor-informed matching. A variant must pass all
filters to receive the PASS designation.

---

## Multi-Reference Support

Multiple reference files can be supplied for joint filtering:

```
kmerdet filter -i results.tsv -r panel_v1.tsv -r panel_v2.tsv -o filtered.tsv
```

References are merged by taking the union of all expected variants. A detected variant
passes if it matches any reference. This supports workflows where the tumor mutation
list is updated over time (e.g., additional mutations discovered in a second biopsy).

---

## Panel-of-Normals Integration

When a PoN database is supplied, an additional filtering layer removes variants
recurrently observed in healthy controls:

```
kmerdet filter -i results.tsv -r reference.tsv --pon normals.pon -o filtered.tsv
```

A variant is filtered as `PON_ARTIFACT` if it appears in more than `--pon-max-freq`
(default: 5%) of normal samples. The PoN catches systematic sequencing errors (GGC
motifs), reference genome errors, common germline variants, and recurrent CHIP mutations.
The false-positive-analysis.md research document details the construction and expected
effectiveness of a PoN.

### PoN Construction

The PoN is built from detection results on a cohort of healthy control samples processed
through the identical pipeline:

```
kmerdet detect -d normal_01.jf -t targets/ -o normals/normal_01.tsv
kmerdet detect -d normal_02.jf -t targets/ -o normals/normal_02.tsv
...
kmerdet pon-build --input normals/ --output panel.pon
```

The PoN database stores `(target, normalized_variant) -> frequency` mappings in a
compact binary format. Storage is minimal: approximately 50 bytes per target-variant
pair.

---

## Streaming Design

The filter operates as a streaming processor: it reads one record at a time from the
input, applies all filter criteria, and writes passing records to the output. It does not
need to hold all results in memory simultaneously.

The reference file and PoN database are loaded into memory at startup (both are small:
a 50-target reference is <10 KB, a PoN is <50 KB). Input records are processed in a
single pass.

This streaming design enables piping:

```
kmerdet detect -d sample.jf -t targets/ | kmerdet filter -r reference.tsv -o filtered.tsv
```

---

## Research Backing

### INDEL Normalization

The Tan et al. (2015) paper "Unified representation of genetic variants" established
the normalization algorithm used by bcftools norm, vt normalize, and GATK
LeftAlignAndTrimVariants. The indel-filtering.md research document shows that applying
this normalization before filtering resolves the representation mismatch problem
identified in the thesis.

### Adaptive Filtering

The adaptive-filtering.md research document proposes coverage-aware thresholds that
replace the fixed `min_vaf`, `min_coverage`, and `min_expression` values. While the
initial implementation uses fixed thresholds, the filter architecture supports plugging
in adaptive thresholds computed from coverage statistics (see the `coverage` subcommand).

### Panel-of-Normals

The PoN approach follows GATK Mutect2's proven pattern. The false-positive-analysis.md
document catalogs the artifact sources a PoN addresses: systematic sequencing errors,
reference genome errors, common germline variants, and CHIP mutations. The thesis
achieved zero false positives without a PoN, but this was on a small cohort with
tumor-informed targets. Scaling to larger panels or de novo discovery will require
PoN filtering.

### CHIP Filtering

CHIP variants are documented as a major confounding factor in cfDNA analysis. Razavi
et al. (Nature Medicine, 2019) found CHIP mutations in 10-30% of cancer patients. The
tumor-informed approach provides inherent protection (CHIP variants are unlikely to match
tumor-specific mutations), but the PoN adds an additional safety layer for any CHIP
variants that coincidentally resemble tumor mutations.

---

## Acceptance Criteria

### Unit Tests

- [ ] Reference mode: exact match on CHROM/POS/REF/ALT/TYPE
- [ ] Reference mode: normalization resolves equivalent INDEL representations
- [ ] Reference mode: chr-prefix normalization (chr17 matches 17)
- [ ] Use-alt mode: exact sequence match
- [ ] Use-alt mode: substring match (detected contains expected)
- [ ] Use-alt mode: overlap match with configurable minimum
- [ ] Normalization: left-alignment of deletion in homopolymer (4 equivalent forms -> 1)
- [ ] Normalization: left-alignment of insertion in repeat
- [ ] Normalization: SNVs pass through unchanged
- [ ] Normalization: complex indels that cannot be simplified remain unchanged
- [ ] Threshold filters: min_vaf, min_expression, min_coverage correctly applied
- [ ] Type exclusion: --exclude-type Reference removes Reference calls
- [ ] Multi-reference: union of two reference files produces correct matches
- [ ] PoN filtering: variants above pon-max-freq are filtered
- [ ] Empty input produces empty output (not an error)
- [ ] Streaming: output appears incrementally, not buffered until completion

### Integration Tests

- [ ] Normalization output matches `bcftools norm` for a test set of INDELs
- [ ] Filter pipeline: detect | filter produces correct filtered output
- [ ] Thesis regression: specific INDEL matching failures from the validation are resolved
- [ ] All output formats preserve filter annotations correctly

### Edge Cases

- [ ] Variant at position 1 (cannot left-align further)
- [ ] Variant at target boundary (limited flanking context)
- [ ] Multiple variants at the same position (multi-allelic)
- [ ] Reference file with duplicate entries (handled gracefully)
- [ ] Detected variant type differs from reference type (e.g., Complex vs Deletion)
