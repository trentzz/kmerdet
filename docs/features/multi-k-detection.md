# Feature: Multi-k Detection -- Multiple K-mer Length Support

## What It Does

Multi-k detection runs variant detection at multiple k-mer lengths (e.g.,
k=21, k=31, k=41) and combines results to improve sensitivity and confidence.
Each k value provides a different view of the data: shorter k-mers catch large
INDELs and low-VAF variants, while longer k-mers provide higher specificity in
repetitive regions. By combining evidence across k values, kmerdet detects
variants that any single k would miss and assigns higher confidence to variants
corroborated at multiple lengths.

## Why It Matters

### The INDEL Sensitivity Crisis

The thesis validation study revealed the single most critical limitation of the
k-mer approach: INDEL sensitivity collapses as INDEL size increases. At k=31:

| INDEL Size | Sensitivity | Junction-spanning k-mers |
|------------|-------------|--------------------------|
| 1-3 bp     | Moderate    | 28-30 of 31              |
| 4-7 bp     | Low         | 24-27 of 31              |
| 10 bp      | Very low    | 21 of 31 (68%)           |
| 15 bp      | ~10%        | 16 of 31 (52%)           |
| 20 bp      | ~3%         | 11 of 31 (35%)           |
| 25+ bp     | Near zero   | 6 or fewer               |

Many clinically actionable cancer mutations are large INDELs: NPM1 4bp insertion
(AML), FLT3-ITD (variable length, often >20 bp), BRCA1/2 frameshifts, EGFR exon
19 deletions (15-24 bp). A tool that achieves only 38% overall INDEL sensitivity
cannot serve as a standalone clinical assay.

At k=21, the same insertions have substantially more junction-spanning k-mers:
a 15 bp insertion retains 6 of 21 junction k-mers (29%), compared to 16 of 31
(52%) at k=31. While k=21 has lower per-k-mer count, the increased junction
coverage compensates at moderate VAF.

### SNV Sensitivity Improvement

The thesis achieved 77% SNV sensitivity at k=31. Multi-k detection can recover
some of the 23% missed SNVs through two mechanisms:

1. **Low-VAF recovery**: At very low VAF (0.05-0.1%), variant k-mer counts may
   fall below the threshold at k=31 but remain above threshold at k=21 (due to
   higher per-k-mer coverage from shorter k-mers in short cfDNA fragments).

2. **Repetitive context**: Some SNVs in near-repetitive regions produce
   ambiguous walks at k=31 but clean walks at k=41, where the longer k-mers
   disambiguate the repeat.

Conservative estimates from the research suggest 5-10% SNV sensitivity
improvement and 15-40% INDEL sensitivity improvement from multi-k detection.

## Research Backing

### K-mer Length Tradeoffs (kmer-length-tradeoffs.md)

The tradeoff analysis establishes that no single k is optimal for all variant
types and coverage levels. At 0.1% VAF with 5000x coverage:

- k=21: expected variant k-mer count = 4.33 per k-mer (above threshold)
- k=31: expected count = 4.00 (above threshold)
- k=43: expected count = 3.60 (marginal at threshold=2)

At 0.05% VAF, k=43 drops below threshold (1.80) while k=21 remains viable
(2.17). The cfDNA fragment size interaction further favors shorter k: tumor-
derived cfDNA fragments are enriched at 90-150 bp, where k=43 wastes 50% of
bases in edge zones versus 24% for k=21.

### Multi-k Strategy (multi-k-strategy.md)

Four approaches were evaluated:

1. **Independent multi-k with union merge**: Run pipeline at each k, take union.
   Simplest implementation, highest resource cost.

2. **Consensus voting**: Weight detection across k values by variant type.
   SNVs weight k=31 highest; long INDELs weight k=21 highest.

3. **Adaptive k per target**: Precompute optimal k per target based on variant
   type and reference characteristics. Walking cost identical to single-k.

4. **Bayesian evidence combination**: Formulate detection as Bayesian inference
   with evidence from multiple k values. Most principled but requires training
   data.

### Assembly Literature Precedent

SPAdes, MEGAHIT, and SKESA all use multi-k strategies in genome assembly.
SPAdes iterates through k=21,33,55,77, using small k for low-coverage regions
and large k for repeat resolution. The analogy to variant detection is direct:
low-VAF variants are analogous to low-coverage regions, and repetitive target
contexts are analogous to assembly repeats.

The BAYSIC framework (Cantarel et al., 2014) demonstrated that Bayesian
combination of calls from multiple variant callers improves both sensitivity
and specificity. Multi-k detection applies the same principle, treating each k
as a different "caller."

## Design Considerations

### Approach: Phased Implementation

**Phase 1 -- Adaptive k per target** (recommended first):

Each target in the panel is assigned an optimal k value based on its
characteristics. This is precomputed once per panel design:

| Target Characteristic | Assigned k | Rationale |
|----------------------|------------|-----------|
| SNV in non-repetitive region | 31 | Standard, well-validated |
| SNV in low-complexity region | 41 | More specificity needed |
| Short INDEL (1-10 bp) | 31 | Standard k works well |
| Medium INDEL (11-20 bp) | 25 | Better junction coverage |
| Long INDEL (>20 bp) | 21 | Maximum junction spanning |
| ITD (any length) | 21 | Tandem duplications need short k |
| Segmental duplication context | 43 | Paralogous region disambiguation |

The target catalog includes a `recommended_k` field:

```toml
[[targets]]
name = "FLT3_ITD"
variant_type = "itd"
recommended_k = 21

[[targets]]
name = "TP53_R175H"
variant_type = "snv"
recommended_k = 31
```

Walking cost is identical to single-k (each target queries one database).
Multiple jellyfish databases are needed (one per distinct k value used), but
this is only a counting-phase cost.

**Phase 2 -- Independent multi-k with consensus voting**:

Run walking at all k values for all targets, then merge results with weighted
consensus voting. The voting weights are variant-type-specific:

For SNVs: w(k=21) = 0.7, w(k=31) = 1.0, w(k=41) = 0.9
For short INDELs: w(k=21) = 1.0, w(k=31) = 0.9, w(k=41) = 0.7
For long INDELs: w(k=21) = 1.0, w(k=31) = 0.5, w(k=41) = 0.2

Consensus score = sum(w(k) * detected(k)) / sum(w(k))

Decision rules:
- consensus >= 0.7: call variant (high confidence)
- consensus >= 0.4: call variant (moderate confidence, flag for review)
- single-k detection with rVAF > 0.5%: call variant (low confidence)
- otherwise: no call

**Phase 3 -- Bayesian evidence combination**:

Formulate detection as posterior inference:

```
P(variant | E_k1, E_k2, E_k3) ~
    P(E_k1 | variant) * P(E_k2 | variant) * P(E_k3 | variant) * P(variant)
```

Requires calibrated detection sensitivity models per variant type per k, fitted
from training data. Phase 2 generates the multi-k evidence needed for training.

### rVAF Reconciliation

Different k values produce different rVAF estimates for the same variant because
the NNLS decomposition operates on different k-mer sets. Three reconciliation
strategies:

1. **Coverage-weighted average**: Weight by median k-mer count at each k, since
   higher coverage produces more accurate rVAF estimates.

2. **Minimum variance**: Use rVAF from the k value with the lowest coefficient
   of variation in k-mer counts (most uniform coverage).

3. **Weighted by variant type**: For SNVs, prefer k=31 rVAF (most validated).
   For INDELs, prefer the k value with the most junction-spanning k-mers.

### Resource Cost

| Resource | Single-k (k=31) | Triple-k (k=21,31,41) | Overhead |
|----------|-----------------|------------------------|----------|
| Counting time | ~1.5 min | ~4.5 min serial, ~1.5 min parallel | 1-3x |
| Walking time | ~2 min | ~6 min serial, ~2 min parallel | 1-3x |
| Memory peak | ~720 MB | ~2.1 GB all loaded, ~768 MB sequential | 1-3x |
| Disk (.jf files) | ~400 MB | ~1.2 GB | 3x |
| Total pipeline | ~6 min | ~7 min parallel, ~12 min serial | 1.2-2x |

For targeted panels (50 targets), multi-k is feasible on standard workstations
with 8-16 GB RAM. Walking can be parallelized across k values per target using
rayon.

### CLI Interface

```bash
# Adaptive k per target (Phase 1) -- targets specify their own k
kmerdet detect -d sample_k21.jf -d sample_k31.jf -t targets/ -o results.tsv

# Multi-k with explicit databases (Phase 2)
kmerdet detect \
    --db sample_k21.jf \
    --db sample_k31.jf \
    --db sample_k41.jf \
    --multi-k-mode consensus \
    -t targets/ -o results.tsv

# Config file approach
# [detect]
# databases = ["sample_k21.jf", "sample_k31.jf", "sample_k41.jf"]
# multi_k_mode = "consensus"
```

### Variant Matching Across k Values

Two calls at different k values match if they describe the same genomic variant.
Matching requires left-alignment normalization since different k values may
produce slightly different breakpoint coordinates for INDELs. The matching
criteria:

- Same target
- Same variant type (SNV, insertion, deletion, ITD)
- Same chromosome and position (after normalization)
- Same ref/alt alleles (after normalization)

Variants that match across k values are merged into a single call with
multi-k annotations. Variants detected at only one k value are retained
with lower confidence.

## Acceptance Criteria

### Phase 1 (Adaptive k per target)

1. Targets with `recommended_k` field are walked using the specified k value
   and the corresponding jellyfish database.

2. Targets without `recommended_k` default to k=31.

3. Multiple jellyfish databases can be loaded simultaneously (one per k value).

4. Detection output includes a `KMER_K` field indicating which k was used.

### Phase 2 (Consensus voting)

5. All targets are walked at all available k values (parallelized).

6. Variant calls are matched across k values using normalized coordinates.

7. Consensus scores are computed with variant-type-specific weights.

8. Detection pattern (which k values detected the variant) is reported.

9. Long INDELs (>15 bp) are detected that were missed at k=31 alone.

10. Specificity is maintained or improved versus single-k detection
    (consensus down-weights single-k false positives).

### Performance

11. Multi-k detection at 3 k values completes in under 2x the wall-clock
    time of single-k detection when parallelized.

12. Memory usage with 3 databases loaded sequentially does not exceed
    1 GB for a typical 50-target panel.

### Validation

13. On the thesis validation dataset, multi-k detection improves SNV
    sensitivity by at least 5 percentage points over single-k (k=31).

14. On the thesis validation dataset, multi-k detection improves INDEL
    sensitivity by at least 15 percentage points over single-k.
