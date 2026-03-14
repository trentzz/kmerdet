# Gene Fusion Detection via K-mer Walking

## Overview

Gene fusions --- chimeric genes formed by the juxtaposition of sequences from two
previously separate genes --- are among the most clinically important structural
variants in oncology. They serve as diagnostic markers, therapeutic targets, and MRD
monitoring anchors. In the context of k-mer-based variant detection, gene fusions
create a distinctive signature: novel junction k-mers at the fusion breakpoint that
are absent from either partner gene's reference sequence.

This document describes how the kmerdet algorithm can detect gene fusions, the target
design strategy required, the implementation plan, and the clinical fusion targets
most relevant to liquid biopsy panels.

## K-mer Signatures of Gene Fusions

### Breakpoint Junction K-mers

A gene fusion joining gene A (at position pA) to gene B (at position pB) creates a
chimeric sequence:

```
Gene A:    ...ACGTACGTACGT|
Gene B:                   |TGCATGCATGCA...
Fusion:    ...ACGTACGTACGT|TGCATGCATGCA...
                          ^
                     Breakpoint
```

The k-mers spanning the breakpoint are novel --- they do not appear in either gene A
or gene B individually. For k=31, there are up to 30 novel junction k-mers (those
containing at least one base from each partner gene). These junction k-mers are the
primary detection signal.

### K-mer Count Patterns

In a sample carrying a fusion at variant allele frequency f:

- **Reference k-mers** (entirely within gene A or gene B): Present at normal coverage
  C. The fusion does not affect these counts.
- **Junction k-mers**: Present at coverage approximately f * C. At 1% VAF with 1000x
  coverage, junction k-mers have ~10 counts.
- **K-mers at the breakpoint boundary** (last k-1 k-mers of gene A before breakpoint):
  Coverage is slightly reduced because the fusion allele contributes reads that
  transition to gene B sequence instead of continuing in gene A.

### How Walking Discovers Fusions

Given a fusion-specific target FASTA containing the chimeric sequence:

1. The walker starts from reference k-mers of gene A (upstream of breakpoint).
2. At the breakpoint, the reference path transitions to gene B sequence.
3. If the sample carries the fusion, junction k-mers have non-zero counts.
4. The walker extends through the junction into gene B k-mers.
5. The graph contains a path through junction k-mers (alternative path) alongside
   the expected reference path.
6. The classifier identifies the event based on the path divergence at the breakpoint.

If the sample does NOT carry the fusion, junction k-mers have zero counts, the walker
cannot extend through the breakpoint, and no alternative path is found. The target
reports as reference-only.

## Target Design for Fusion Detection

### Basic Chimeric Target Construction

Each fusion target is a FASTA sequence representing the expected chimeric transcript
or genomic sequence:

```fasta
>BCR_ABL1_e13a2
ATGTTGGAGATCTGCCTGAAGCTGGTGGGCTGCAAATCCGTACTGACATCCGTGGAGCTG
CATTCATGTGACTTTGATGGTGACTTTGTCTTCCTCATCTCCCCTCTCTTCCTCATTTTTC
ACGGATCACCATATTCCTTCCAAGGTCATTGATGACTCTACTCAGAAACTCAATGATGATG
ATGATGATGATGATGATGATGATCTACATTGATAATGAGATGGAGTACCTGAGACGGGGCA
```

The target contains:
- Upstream sequence from gene A (at least k-1 bases, ideally 100-200bp for context)
- Downstream sequence from gene B (at least k-1 bases, ideally 100-200bp)
- The junction point where the two sequences meet

### Enumerating Breakpoint Positions

For a given gene fusion pair, multiple breakpoint positions may be biologically valid:

```
Exon-exon fusions (most common in expressed fusions):
  Gene A exon 13 | Gene B exon 2    (e13a2, "major" BCR-ABL)
  Gene A exon 14 | Gene B exon 2    (e14a2, "major" BCR-ABL)
  Gene A exon 1  | Gene B exon 2    (e1a2, "minor" BCR-ABL)
```

Each breakpoint position requires a separate target FASTA. For a fusion with N
possible breakpoint positions, N targets are needed.

### Target Naming Convention

A consistent naming convention aids downstream analysis:

```
>{GENE_A}_{GENE_B}_{breakpoint_id}
```

Examples:
```
>BCR_ABL1_e13a2
>BCR_ABL1_e14a2
>BCR_ABL1_e1a2
>EML4_ALK_v1_e13a20
>EML4_ALK_v2_e20a20
>EML4_ALK_v3a_e6a20
```

### Target Length Considerations

The target should be long enough to:
- Provide sufficient reference k-mers for walking context (at least 100bp on each side
  of the breakpoint)
- Allow the graph to fully resolve the junction region
- Avoid edge effects from the target boundaries

Recommended target length: 200-400bp total (100-200bp from each partner gene).

For liquid biopsy where cfDNA fragments are ~170bp, the target need not be longer
than ~300bp, as no single fragment will provide more than ~170bp of contiguous
sequence.

### Handling Alternative Splicing

Alternative splicing in one or both partner genes can create multiple valid fusion
transcripts:

```
EML4-ALK variant 1: EML4 exon 13 | ALK exon 20
EML4-ALK variant 2: EML4 exon 20 | ALK exon 20
EML4-ALK variant 3a: EML4 exon 6  | ALK exon 20  (partial exon 6)
EML4-ALK variant 3b: EML4 exon 6  | ALK exon 20  (full exon 6)
```

Each variant produces different junction k-mers and requires a separate target. A
comprehensive fusion panel includes all known clinically relevant splice variants.

### Reciprocal Fusion Targets

Some fusions produce detectable reciprocal products:

```
Forward:    BCR exon 13 | ABL1 exon 2   (BCR-ABL1)
Reciprocal: ABL1 exon 1 | BCR exon 14   (ABL1-BCR)
```

In most clinical contexts, only one orientation is therapeutically relevant, but
detecting the reciprocal fusion can increase sensitivity and confirm the event. Both
orientations should be included in a comprehensive panel.

## Challenges in Fusion Detection

### Multiple Possible Breakpoint Positions

For well-characterized fusions, the number of known breakpoint positions is manageable
(2-5 per fusion pair). However, for less well-characterized fusions, the breakpoint
may occur at any position within an exon, requiring many targets or a sliding-window
approach.

**Mitigation:** Focus on canonical exon-exon boundaries for initial implementation.
Allow a "fuzzy" mode that generates targets at regular intervals across exon
boundaries for discovery.

### Cryptic Breakpoints

Some fusions have breakpoints in intronic or intergenic regions, far from canonical
exon boundaries. These create fusion transcripts that include partial intronic
sequence.

**Mitigation:** For known recurrent fusions with non-canonical breakpoints (e.g., some
FGFR2 fusions), include these specific breakpoint positions in the target panel. For
novel breakpoints, a discovery mode using de novo assembly of unmapped k-mers would be
needed (out of scope for initial implementation).

### Repetitive Sequences Near Breakpoints

If the sequences flanking a breakpoint are repetitive (e.g., Alu elements, LINE
elements), the junction k-mers may not be unique in the genome, leading to false
positive counts from other genomic loci.

**Mitigation:**
- Use longer k values (k=31 or k=51) to increase specificity of junction k-mers.
- Check junction k-mer uniqueness against the reference genome during target design.
- Flag targets where junction k-mers have high reference genome frequency.

### Low VAF Detection

At the VAF levels relevant to MRD monitoring (0.01% to 1%), junction k-mers may
have very low counts:

| Coverage | VAF   | Expected junction k-mer count |
|----------|-------|-------------------------------|
| 1000x    | 1%    | ~10                           |
| 1000x    | 0.1%  | ~1                            |
| 1000x    | 0.01% | ~0.1 (undetectable)           |
| 10000x   | 0.1%  | ~10                           |
| 10000x   | 0.01% | ~1                            |

At very low VAF, junction k-mers are indistinguishable from sequencing error k-mers
by count alone. The walking threshold (default: count >= 2, ratio >= 0.05) filters
many true positives at low VAF.

**Mitigation:**
- Lower walking thresholds for fusion-specific targets (e.g., count >= 1).
- Use UMI-deduplicated counts where each count represents an independent molecule.
- Require multiple junction k-mers to support the fusion call (consensus across
  the junction region).
- Apply fusion-specific statistical models (binomial probability of observing N
  junction k-mers by chance).

### Novel Fusion Partners

Fusions involving genes not in the target panel cannot be detected. This is a
fundamental limitation of the targeted k-mer approach.

**Mitigation:**
- Comprehensive panels covering all recurrent fusions for the relevant cancer type.
- For discovery, complement k-mer analysis with other methods (RNA-seq fusion callers,
  whole-genome SV callers).

## Common Fusion Targets for Liquid Biopsy Panels

### Hematologic Malignancies

| Fusion           | Disease      | Frequency | Breakpoint Variants | Notes                            |
|------------------|-------------|-----------|---------------------|----------------------------------|
| BCR-ABL1         | CML, ALL    | ~95% CML  | e13a2, e14a2, e1a2  | Defining lesion for CML          |
| PML-RARA         | APL (AML)   | ~98% APL  | bcr1, bcr2, bcr3    | Critical for APL diagnosis       |
| CBFB-MYH11       | AML         | ~8% AML   | Type A (most common) | Core-binding factor AML          |
| RUNX1-RUNX1T1    | AML         | ~7% AML   | Standard             | Core-binding factor AML          |
| KMT2A-MLLT3      | AML, ALL    | ~5% AML   | Multiple             | MLL-rearranged leukemia          |
| NPM1-ALK         | ALCL        | ~75% ALCL | Standard             | Anaplastic large cell lymphoma   |
| IGH-BCL2         | FL          | ~85% FL   | MBR, mcr             | Follicular lymphoma              |
| EWSR1-FLI1       | Ewing       | ~85%      | Type 1, Type 2       | Ewing sarcoma                    |

### Solid Tumor Fusions

| Fusion           | Disease      | Frequency    | Breakpoint Variants   | Notes                        |
|------------------|-------------|-------------|-----------------------|------------------------------|
| EML4-ALK         | NSCLC       | ~3-5%       | v1, v2, v3a/b         | ALK inhibitor responsive     |
| ROS1 fusions     | NSCLC       | ~1-2%       | Multiple partners     | Crizotinib responsive        |
| RET fusions      | NSCLC, MTC  | ~1-2%       | KIF5B-RET, CCDC6-RET  | Selpercatinib responsive     |
| NTRK fusions     | Multiple    | <1% (most)  | Many partners         | Larotrectinib, entrectinib   |
| FGFR2 fusions    | Cholang.    | ~10-15%     | Multiple partners     | Pemigatinib responsive       |
| FGFR3 fusions    | Bladder     | ~3%         | FGFR3-TACC3          | Erdafitinib responsive       |
| TMPRSS2-ERG      | Prostate    | ~50%        | Multiple              | Diagnostic, not therapeutic  |
| PAX3/7-FOXO1     | RMS         | ~80%        | PAX3, PAX7 variants  | Rhabdomyosarcoma             |

## Implementation Plan for kmerdet

### Phase 1: Fusion Target FASTA Support (Already Functional)

The current kmerdet `detect` subcommand already supports fusion detection when given
appropriate target FASTAs. No additional code changes are needed for basic fusion
detection.

**Verification steps:**
1. Create a test fusion target FASTA (e.g., BCR-ABL1 e13a2).
2. Create a mock k-mer database with junction k-mers at known counts.
3. Run `kmerdet detect` and verify that the fusion is detected and classified.
4. Verify that rVAF is correctly computed for the fusion allele.

### Phase 2: Fusion Target Generator Tool

A new subcommand or utility to automatically generate fusion target FASTAs from
gene pair specifications and a reference genome.

**Input:**
```toml
[fusion_targets]
reference_genome = "hg38.fa"
reference_gtf = "gencode.v38.annotation.gtf"

[[fusion_targets.fusions]]
gene_a = "BCR"
gene_b = "ABL1"
breakpoints = ["e13a2", "e14a2", "e1a2"]
flank_size = 150

[[fusion_targets.fusions]]
gene_a = "EML4"
gene_b = "ALK"
breakpoints = ["v1_e13a20", "v2_e20a20", "v3a_e6a20"]
flank_size = 150
```

**Output:** Directory of FASTA files, one per fusion breakpoint variant.

**Implementation approach:**
1. Parse the reference GTF to extract exon coordinates for each gene.
2. Extract exon sequences from the reference genome FASTA.
3. For each specified breakpoint, concatenate the appropriate exon sequences.
4. Write the chimeric sequence as a FASTA file with descriptive headers.
5. Validate that junction k-mers are unique in the reference genome (flag
   non-unique junctions).

**Rust modules involved:**
- New `src/fusion/` module for fusion target generation.
- Uses `needletail` for FASTA reading (already a dependency).
- GTF parsing via a lightweight parser (e.g., `noodles-gtf` or custom).

### Phase 3: Fusion-Specific Reporting

Enhanced output for fusion variant calls:

```tsv
target_name        variant_type  gene_a  gene_b  breakpoint_a  breakpoint_b  rVAF    confidence
BCR_ABL1_e13a2     FUSION        BCR     ABL1    exon13:end    exon2:start   0.015   35.2
EML4_ALK_v1        FUSION        EML4    ALK     exon13:end    exon20:start  0.008   28.7
```

**Fields added for fusion calls:**
- `gene_a`, `gene_b`: Partner gene names (parsed from target name).
- `breakpoint_a`, `breakpoint_b`: Breakpoint positions within each gene.
- `fusion_variant`: Canonical fusion variant name (e.g., "e13a2").
- `junction_kmer_count`: Number of junction k-mers detected.
- `junction_kmer_fraction`: Fraction of expected junction k-mers that were observed.

**VCF representation:**
```vcf
chr9  5073770  BCR_ABL1_e13a2  N  <FUS>  35.2  PASS  SVTYPE=BND;MATEID=ABL1_BCR;GENE_A=BCR;GENE_B=ABL1;FUSION_VAR=e13a2  GT:VAF  0/1:0.015
```

### Phase 4: Fusion-Specific Confidence Scoring

Fusion calls benefit from a specialized confidence model:

- **Junction k-mer consensus:** The fraction of expected junction k-mers that are
  observed. A true fusion should produce most of the expected ~k-1 junction k-mers.
  Low consensus suggests a false positive from a single error k-mer.
- **Strand consistency:** For strand-specific libraries, junction k-mers should appear
  on the expected strand. Inconsistent strand suggests artifact.
- **Partner gene balance:** Both gene A and gene B reference k-mers should be at
  normal coverage. Abnormal coverage in one partner may indicate a different event
  (e.g., gene deletion, not fusion).
- **Binomial model:** Given reference coverage C, the probability of observing N
  junction k-mers by chance from sequencing error can be computed as:
  ```
  P(N | error) = Binom(N; K, epsilon^k)
  ```
  where K is the number of expected junction k-mers and epsilon is the per-base
  error rate.

## VAF Interpretation for Fusions

The rVAF computed by the NNLS quantification represents the fusion allele frequency:

- **rVAF = fusion reads / total reads at the locus**
- For a sample with 1% tumor fraction carrying a heterozygous fusion:
  rVAF ~ 0.5% (fusion is on one allele of the tumor cells).
- For a sample with 1% tumor fraction carrying a homozygous fusion (or fusion in
  CML where BCR-ABL is the driver): rVAF ~ 1%.

In practice, fusion VAF in liquid biopsy correlates with tumor burden and is used
for MRD monitoring:

- **Diagnosis:** High fusion VAF confirms the fusion and disease.
- **Treatment response:** Decreasing fusion VAF indicates response.
- **MRD monitoring:** Persistent low-level fusion VAF indicates residual disease.
- **Relapse detection:** Rising fusion VAF from undetectable indicates relapse.

The sensitivity limit for fusion detection is determined by:
- Sequencing depth (more reads = more junction k-mers)
- UMI deduplication efficiency (each UMI family represents one original molecule)
- K-mer counting threshold (minimum count to call a k-mer as present)

With 10,000x deduplicated coverage and a minimum count of 2, the theoretical
detection limit is approximately 0.02% VAF (2 independent molecules out of 10,000).

## References

- Mitelman, F., Johansson, B., & Mertens, F. (2007). The impact of translocations
  and gene fusions on cancer causation. Nature Reviews Cancer, 7(4), 233-245.
- Parker, B.C., & Zhang, W. (2013). Fusion genes in solid tumors: an emerging target
  for cancer diagnosis and treatment. Chinese Journal of Cancer, 32(11), 594-603.
- Heyer, E.E., et al. (2019). Diagnosis of fusion genes using targeted RNA sequencing.
  Nature Communications, 10(1), 1388.
- Soda, M., et al. (2007). Identification of the transforming EML4-ALK fusion gene
  in non-small-cell lung cancer. Nature, 448(7153), 561-566.
- Druker, B.J., et al. (2001). Efficacy and safety of a specific inhibitor of the
  BCR-ABL tyrosine kinase in chronic myeloid leukemia. New England Journal of
  Medicine, 344(14), 1031-1037.
- Cocco, E., Scaltriti, M., & Drilon, A. (2018). NTRK fusion-positive cancers and
  TRK inhibitor therapy. Nature Reviews Clinical Oncology, 15(12), 731-747.
- Abou-Elkacem, L., et al. (2023). Liquid biopsy for fusion gene detection:
  current status and future directions. Clinical Chemistry, 69(6), 565-576.
