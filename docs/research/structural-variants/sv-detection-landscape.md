# Structural Variant Detection via K-mer Methods

## Overview

Structural variants (SVs) represent genomic rearrangements typically larger than 50bp,
including deletions, insertions, duplications, inversions, translocations, and gene
fusions. In the liquid biopsy setting, circulating tumor DNA (ctDNA) carries these
rearrangements as fragmented cell-free DNA (cfDNA), typically 150-170bp in length.
Detecting SVs from k-mer databases presents unique challenges and opportunities
compared to SNV and small INDEL detection.

The kmerdet algorithm --- load target FASTAs, walk k-mers via DFS, build directed
weighted graph, find alternative paths, classify variants, quantify with NNLS --- was
designed primarily for SNV and small INDEL detection. This document surveys how each
class of structural variant manifests in k-mer space, assesses the current capability
of the km/kmerdet approach, and identifies the key challenges for extending detection
to larger and more complex rearrangements.

## SV Types Relevant to Liquid Biopsy

### Large Deletions (>50bp)

**K-mer space manifestation:**
- K-mers spanning the deleted region are absent from the sample database (count = 0).
- K-mers flanking the deletion are present at normal coverage.
- A novel "junction k-mer" exists at the breakpoint, composed of sequence from the
  left flank and right flank of the deletion. This junction k-mer is not present in
  the reference.
- For a deletion of size D with k-mer length k:
  - Exactly D reference k-mers are missing (those entirely within the deleted region).
  - Up to k-1 reference k-mers have reduced counts (those partially overlapping the
    breakpoint).
  - Up to k-1 novel junction k-mers are created at the breakpoint.

**Current capability:**
- Deletions smaller than k-1 bp (30bp at k=31) are well-handled by the standard
  walking algorithm. The walker traverses from reference k-mers, discovers the
  alternative (shorter) path through the graph, and the classifier identifies the
  deletion.
- Deletions between k-1 and read_length - k are partially detectable: junction k-mers
  exist in the database (reads span the breakpoint), but the graph may become
  disconnected because no single k-mer bridges the gap in reference space.
- Deletions larger than read_length are undetectable from k-mer counts alone, as no
  sequencing read spans the breakpoint.

**Key challenges:**
- Graph disconnection for deletions >k-1 bp. The reference k-mer path breaks, and the
  walker cannot traverse from one side of the deletion to the other through normal
  k-mer overlap.
- Target design must anticipate the deletion. The current approach requires a target
  FASTA that includes the region to be deleted; if the deletion extends beyond the
  target, it is invisible.
- At low VAF, junction k-mers may have counts below the walking threshold, making them
  indistinguishable from sequencing errors.

### Large Insertions (>k bp)

**K-mer space manifestation:**
- K-mers from the inserted sequence are novel (not in reference). For insertions larger
  than k-1 bp, some of these novel k-mers are entirely internal to the insertion and
  have no overlap with reference k-mers.
- Junction k-mers exist at both ends of the insertion, connecting reference sequence to
  inserted sequence.
- For an insertion of size I:
  - 2*(k-1) junction k-mers span the two breakpoints.
  - I - k + 1 novel internal k-mers exist within the insertion.
  - Reference k-mers flanking the insertion site remain at normal coverage.

**Current capability:**
- Insertions smaller than k-1 bp (30bp at k=31) are well-detected. The walker
  discovers the alternative (longer) path through novel k-mers, and the classifier
  identifies the insertion.
- Insertions between k-1 and read_length - k are partially detectable. Junction k-mers
  exist, but internal insertion k-mers create a region the walker cannot reach from
  reference k-mers alone. The graph has a "bridge" of novel k-mers disconnected from
  the main reference path.
- Insertions larger than read_length are undetectable, as no read contains both a
  reference k-mer and an insertion k-mer.

**Key challenges:**
- Disconnected graph regions for large insertions. The novel internal k-mers form an
  island unreachable from the reference-seeded walk.
- Distinguishing true insertion k-mers from sequencing error k-mers is harder when the
  insertion is large (more novel k-mers to validate).
- Insertion size estimation requires traversing the entire inserted sequence, which may
  fail if coverage is uneven across the insertion.

### Internal Tandem Duplications (ITDs)

**K-mer space manifestation:**
- The duplicated region produces k-mers at elevated counts (approximately 2x for
  heterozygous duplication at 100% tumor fraction).
- Novel junction k-mers are created at the duplication boundary, where the end of the
  duplicated region connects back to its beginning.
- For an ITD of size D:
  - D reference k-mers have elevated counts.
  - Up to k-1 novel junction k-mers exist at the tandem junction.

**Current capability:**
- ITDs up to k-1 bp are well-modeled. The walker discovers the loop path through the
  graph, and the existing ITD classifier handles this case.
- FLT3-ITD is a primary use case documented in the thesis. FLT3-ITDs range from 3 to
  >400bp, with a median of ~30bp. At k=31, ITDs up to 30bp are within the standard
  detection capability.
- The thesis reports successful FLT3-ITD detection as a key validation case.

**Key challenges:**
- Large ITDs (>k-1 bp) create graph structures similar to large insertions, with
  disconnected novel junction k-mers.
- ITD size heterogeneity: a single sample may harbor multiple ITD alleles of different
  sizes, creating a complex graph with multiple alternative paths.
- At low VAF, the elevated k-mer counts from the duplication may not be distinguishable
  from normal coverage variation.

### Gene Fusions

**K-mer space manifestation:**
- Novel junction k-mers at the fusion breakpoint, composed of sequence from two
  different genes (potentially on different chromosomes).
- K-mers from both partner genes are present at expected coverage.
- For a fusion with breakpoint at position B:
  - Up to k-1 novel junction k-mers span the fusion point.
  - K-mers from gene A upstream of the breakpoint are present.
  - K-mers from gene B downstream of the breakpoint are present.

**Current capability:**
- Gene fusions are detectable with the current km/kmerdet approach, provided that
  appropriate fusion-specific target FASTAs are supplied.
- The target FASTA must contain the chimeric sequence: upstream gene exon concatenated
  with downstream gene exon.
- km already supports fusion detection with this target design strategy.
- The walker discovers the junction k-mers, the graph shows the alternative path
  through the fusion junction, and the classifier can identify the event.

**Key challenges:**
- Target design burden: Every possible fusion must be pre-specified. Novel fusions
  outside the target panel are invisible.
- Multiple possible breakpoint positions within exons require multiple target
  sequences per fusion pair.
- Alternative splicing creates multiple valid fusion transcripts.
- Reciprocal fusions (A-B and B-A) need separate targets.
- Cryptic breakpoints outside canonical exon-exon boundaries are missed by
  exon-boundary-based target design.

### Translocations

**K-mer space manifestation:**
- Essentially identical to gene fusions from a k-mer perspective. The breakpoint
  creates novel junction k-mers joining sequences from two genomic loci.
- Balanced translocations produce two sets of junction k-mers (one for each derivative
  chromosome).

**Current capability:**
- Same as gene fusions: detectable with appropriate target FASTAs.
- No additional algorithmic work is needed beyond what fusion detection requires.

**Key challenges:**
- Same as gene fusions, plus:
  - Balanced translocations require targets for both derivative chromosomes.
  - Breakpoints in repetitive regions create ambiguous junction k-mers.
  - Translocations involving large distances (megabases) are conceptually identical to
    fusions but may involve different biological significance.

### Copy Number Variations (CNVs)

**K-mer space manifestation:**
- CNV gains: K-mers in the amplified region show elevated counts proportional to copy
  number. A region at 4 copies (vs. 2 normal) shows approximately 2x k-mer counts.
- CNV losses: K-mers in the deleted region show reduced counts. Heterozygous deletion
  shows approximately 0.5x counts; homozygous deletion shows near-zero counts.
- No novel junction k-mers are created (unlike deletions with breakpoints).
- CNV boundaries may or may not produce junction k-mers depending on the mechanism
  (tandem duplication vs. dispersed duplication).

**Current capability:**
- Not supported by the current single-target walking approach. The km algorithm is
  designed to find alternative paths in a graph, not to measure regional depth
  changes.
- The `kmerdet coverage` subcommand can report k-mer depth across a target region,
  which is a building block for CNV detection.

**Key challenges:**
- CNV detection requires a fundamentally different approach: comparing k-mer depth
  ratios across regions rather than finding alternative graph paths.
- GC bias, mappability, and library complexity affect k-mer counts and confound
  CNV calling.
- In liquid biopsy, ctDNA fraction modulates the expected fold-change, making CNV
  detection at low tumor fractions extremely difficult.
- Requires a reference (normal) k-mer database for comparison, or a panel-of-normals
  approach.

## Summary Table

| SV Type            | Size Range        | Current Support       | Required Approach                          | Priority for Liquid Biopsy |
|--------------------|-------------------|-----------------------|--------------------------------------------|----------------------------|
| Small deletion     | 1 to k-1 bp      | Full                  | Standard walking + classification          | High (common)              |
| Large deletion     | k to read_len bp  | Partial               | Junction k-mer targets, coverage drop      | Medium                     |
| Very large del.    | >read_len bp      | None                  | Paired k-mer analysis, external SV calls   | Low                        |
| Small insertion    | 1 to k-1 bp      | Full                  | Standard walking + classification          | High (common)              |
| Large insertion    | k to read_len bp  | Partial               | Multi-k, junction k-mer search             | Medium                     |
| ITD (small)        | 3 to k-1 bp      | Full                  | Standard walking + ITD classifier          | High (FLT3-ITD)            |
| ITD (large)        | k to 400bp        | Partial               | Multi-k, dedicated ITD targets             | High (FLT3-ITD)            |
| Gene fusion        | N/A (junction)    | Full (with targets)   | Fusion-specific target FASTAs              | Very High (hematologic)    |
| Translocation      | N/A (junction)    | Full (with targets)   | Same as fusion                             | High (hematologic)         |
| CNV gain           | >1kb              | None                  | K-mer depth ratio analysis                 | Medium (amplifications)    |
| CNV loss           | >1kb              | None                  | K-mer depth ratio analysis                 | Low (redundant with del.)  |

## Prioritization for MRD Monitoring

In the context of minimal residual disease (MRD) monitoring via liquid biopsy, the
priority ranking for SV detection is driven by clinical actionability and prevalence:

### Very High Priority
- **Gene fusions in hematologic malignancies**: BCR-ABL (CML), PML-RARA (APL),
  CBFB-MYH11 (AML), RUNX1-RUNX1T1 (AML). These are defining molecular events used
  for MRD monitoring. Fusion persistence after treatment indicates residual disease.
  The km approach already supports these with appropriate targets.

### High Priority
- **FLT3-ITD**: Most common mutation in AML (~25% of cases). ITD size ranges from
  3 to >400bp. Small ITDs are well-detected; larger ITDs need multi-k or dedicated
  targets.
- **Small INDELs (<30bp)**: Common across many cancer types. Well-supported by current
  approach but sensitivity drops at low VAF.
- **Solid tumor fusions**: EML4-ALK, ROS1, RET, NTRK fusions in NSCLC and other
  solid tumors. Increasingly important for therapy selection and MRD monitoring.

### Medium Priority
- **Large INDELs (30-100bp)**: Less common but clinically significant. Multi-k
  approach is the primary strategy for improvement.
- **Large deletions**: Important for some tumor suppressor losses but harder to detect
  from cfDNA.

### Lower Priority
- **CNVs**: Important for some contexts (HER2 amplification, MYC amplification) but
  require fundamentally different methods and are less suitable for k-mer walking.
- **Very large SVs (>read_length)**: Undetectable from k-mer databases alone;
  require structural information from read pairs or long reads.

## Relationship to Multi-k Strategy

The multi-k approach documented in `docs/research/kmer-length-optimization/` is the
single most impactful improvement for SV detection across multiple categories:

- **k=21**: Extends the "k-1 barrier" from 30bp to 20bp, but captures more INDELs
  in the 20-30bp range that were previously missed. Also reduces the minimum SV size
  that causes graph disconnection.
- **k=31**: Maintains high specificity for SNVs and provides good sensitivity for
  SVs up to 30bp.
- **k=51**: Could improve specificity for fusion detection by requiring longer
  junction matches, reducing false positives from repetitive sequences.

The consensus merging strategy already implemented in kmerdet naturally extends to
SV detection: run detection at multiple k values, merge results, and use agreement
across k values as a confidence signal.

## References

- Audoux, J., et al. (2017). DE-kupl: exhaustive capture of biological variation in
  RNA-seq data through k-mer decomposition. Genome Biology, 18(1), 243.
- Kokot, M., Dlugosz, M., & Deorowicz, S. (2017). KMC 3: counting and manipulating
  k-mer statistics. Bioinformatics, 33(17), 2759-2761.
- Marchet, C., et al. (2020). Data structures based on k-mers for querying large
  collections of sequencing data sets. Genome Research, 31(1), 1-12.
- Nordborg, M., & Tavaré, S. (2002). Linkage disequilibrium: what history has to
  tell us. Trends in Genetics, 18(2), 83-90.
- Layer, R.M., et al. (2014). LUMPY: a probabilistic framework for structural variant
  discovery. Genome Biology, 15(6), R84.
- Chen, X., et al. (2016). Manta: rapid detection of structural variants and indels
  for germline and cancer sequencing applications. Bioinformatics, 32(8), 1220-1222.
- Cameron, D.L., et al. (2017). GRIDSS: sensitive and specific genomic rearrangement
  detection using positional de Bruijn graph assembly. Genome Research, 27(12),
  2050-2060.
- Rausch, T., et al. (2012). DELLY: structural variant discovery by integrated
  paired-end and split-read analysis. Bioinformatics, 28(18), i333-i339.
