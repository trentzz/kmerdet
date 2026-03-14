# Improving Large INDEL Detection (>7bp)

## The Problem

The thesis evaluating the km-based variant detection pipeline reports dramatically
reduced sensitivity for larger INDELs:

| INDEL Size | Sensitivity |
|------------|------------|
| 1-3bp      | ~55%       |
| 4-7bp      | ~38%       |
| >7bp       | ~3%        |

Overall INDEL sensitivity was 38% compared to 77% for SNVs. The >7bp category
is essentially non-functional, representing a critical gap for clinically important
variants such as FLT3-ITDs (median ~30bp), CALR exon 9 INDELs (5-52bp), NPM1 exon
12 insertions (4bp, well-detected), and EGFR exon 19 deletions (9-24bp).

This document analyzes why large INDELs are hard to detect in k-mer space and
presents concrete strategies for improvement.

## Why Large INDELs Are Hard in K-mer Space

### The K-1 Barrier

The fundamental constraint is the k-mer overlap length: two adjacent k-mers share
k-1 bases. For a variant to be fully captured by the walking algorithm, the variant
region must be bridgeable by k-mer overlaps. This creates a hard size limit:

**Maximum INDEL size detectable by single-step k-mer overlap: k-1 bases.**

For k=31, this is 30bp. For k=21, this is 20bp. Here is why:

#### Insertions

Consider an insertion of size I at position P in the reference:

```
Reference: ...AAAA[P]BBBB...
Mutant:    ...AAAA[XXXXXXX]BBBB...    (I bases inserted)
```

K-mers spanning the insertion junction contain bases from both flanking regions (A
and B) and the inserted sequence (X). For a k-mer to "bridge" from the reference
context on one side to the other:

- The k-mer must contain at least 1 base from region A and at least 1 base from
  region B, with up to k-2 bases from the insertion.
- Therefore, the maximum insertion size that allows a bridging k-mer is k-2 bases.
- In practice, the walking algorithm needs a chain of overlapping k-mers through the
  insertion, so the effective limit is that the insertion must be traversable by
  k-mers that share k-1 overlaps with the flanking reference context.

For insertions of size I > k-1:
- Some junction k-mers (at the insertion boundaries) still share overlap with
  reference k-mers.
- But internal insertion k-mers (those entirely within the inserted sequence) form
  an isolated chain with no connection to reference k-mers.
- The walker, starting from reference k-mers, can reach the start of the insertion
  but cannot traverse through it to reach the other side if the insertion is too long.

#### Deletions

Consider a deletion of size D at position P:

```
Reference: ...AAAA[YYYYYYYY]BBBB...   (D bases present)
Mutant:    ...AAAA[--------]BBBB...   (D bases deleted)
           ...AAAABBBB...
```

The deletion creates a junction k-mer composed of flanking sequence:

```
Junction k-mer: ...AAA|BBB...    (k-1 bases from A + rest from B, or vice versa)
```

For deletions of size D > k-1:
- The junction k-mer is entirely novel --- it does not overlap with any reference
  k-mer by k-1 bases, because the reference k-mers adjacent to the deletion include
  bases from the deleted region that are absent in the junction k-mer.
- The walker cannot reach the junction k-mer from reference k-mers on either side
  of the deletion.
- The graph becomes disconnected at the deletion site.

### Graph Disconnection

When the k-1 barrier is exceeded, the k-mer graph exhibits a characteristic
disconnection pattern:

```
For a large deletion (D > k-1):

Reference path:  [A1]-[A2]-...-[Am]---[D1]-[D2]-...-[Dn]---[B1]-[B2]-...-[Bp]
                                       ^deleted region^
Mutant path:     [A1]-[A2]-...-[Am]---[J1]-[J2]-...-[Jq]---[B1]-[B2]-...-[Bp]
                                       ^junction k-mers^

If D > k-1:
  - [Am] and [J1] do NOT share a k-1 overlap
  - [Jq] and [B1] do NOT share a k-1 overlap
  - The walker starting from reference k-mers [Ai] cannot reach [Ji]
  - The walker starting from reference k-mers [Bi] cannot reach [Ji]
  - Junction k-mers [Ji] form an island
```

This disconnection is the root cause of the near-zero sensitivity for large INDELs.

### Coverage Effects at Low VAF

Even for INDELs within the k-1 barrier, low VAF compounds the detection challenge:

- Junction k-mers have counts proportional to VAF * coverage.
- At 0.1% VAF and 1000x coverage, junction k-mers have ~1 count.
- The walking threshold (default: count >= 2, ratio >= 0.05) filters these.
- Larger INDELs affect more junction k-mers, each at low count, making the
  statistical signal weaker per-k-mer but potentially stronger in aggregate.

## Strategy 1: Multi-k Approach (Primary)

### Rationale

Using multiple k-mer lengths shifts the k-1 barrier for each k value:

| k value | k-1 barrier | Detectable INDEL range |
|---------|-------------|----------------------|
| k=21    | 20bp        | 1-20bp               |
| k=31    | 30bp        | 1-30bp               |
| k=41    | 40bp        | 1-40bp               |
| k=51    | 50bp        | 1-50bp               |

By running detection at k=21 and k=31, INDELs up to 30bp are detectable by at least
one k value. Adding k=41 extends this to 40bp, covering the majority of clinically
relevant INDELs.

### Implementation in kmerdet

The multi-k strategy is already partially implemented in kmerdet:

1. **Multiple jellyfish databases:** Create jellyfish databases at each k value from
   the same sequencing data.
   ```bash
   jellyfish count -m 21 -s 1G -C reads.fa -o counts_k21.jf
   jellyfish count -m 31 -s 1G -C reads.fa -o counts_k31.jf
   ```

2. **Parallel detection:** Run `kmerdet detect` independently at each k value.
   ```bash
   kmerdet detect --db counts_k21.jf --targets targets/ -k 21 -o results_k21.tsv
   kmerdet detect --db counts_k31.jf --targets targets/ -k 31 -o results_k31.tsv
   ```

3. **Result merging:** Merge results across k values using consensus logic.
   ```bash
   kmerdet merge --inputs results_k21.tsv results_k31.tsv -o merged.tsv
   ```

### Merging Strategy for Multi-k Results

The merge must handle several cases:

**Case 1: Variant detected at all k values.**
- High confidence. Use the result from the k value with the best rVAF precision
  (typically the larger k, which has fewer false positive branches).
- Consensus across k values is itself a strong confidence signal.

**Case 2: Variant detected at small k only.**
- The variant may be a large INDEL that exceeds the k-1 barrier for larger k values.
- Accept the result from the smaller k, but flag it as "single-k detection."
- Verify that the INDEL size is consistent with the detection pattern (e.g., a 25bp
  deletion detected at k=21 but not k=31 is consistent; a 5bp deletion detected at
  k=21 but not k=31 is suspicious).

**Case 3: Variant detected at large k only.**
- The variant is likely an SNV or small INDEL with better specificity at larger k.
- Accept the result from the larger k.
- This can happen if the smaller k has too many branching paths (reduced specificity)
  that obscure the true variant.

**Case 4: Contradictory results.**
- Different variant types called at different k values for the same target.
- Flag for manual review. May indicate a complex variant.

### Expected Improvement

Based on the k-1 barrier analysis:

| INDEL Size | k=31 only | k=21 + k=31 | k=21 + k=31 + k=41 |
|------------|-----------|-------------|---------------------|
| 1-7bp      | ~38%      | ~55%*       | ~55%                |
| 8-20bp     | ~10%      | ~40%*       | ~40%                |
| 21-30bp    | ~3%       | ~3%         | ~25%*               |
| 31-40bp    | 0%        | 0%          | ~15%*               |

*Estimated. Actual improvement depends on coverage, VAF, and target design.

The improvement at smaller k values comes from two factors:
1. Shorter k-mers are more likely to be covered by reads (each read contributes more
   k-mers at smaller k).
2. The k-1 barrier is lower, allowing the graph to remain connected for larger INDELs.

The tradeoff: smaller k values produce more false positive branches (more k-mers
match by chance), so specificity decreases. The multi-k merge addresses this by
using the larger k for specificity when available.

### Sensitivity vs. Specificity Tradeoff

```
Specificity:  k=51 > k=41 > k=31 > k=21
                (fewer false branches)

Sensitivity:  k=21 > k=31 > k=41 > k=51
(large INDEL)  (larger k-1 barrier)

Optimal strategy: Use the largest k that detects the variant.
```

## Strategy 2: Split-Read K-mer Analysis

### Concept

For very large INDELs (>50bp), individual sequencing reads may span the breakpoint
(split reads). These reads contain junction k-mers that are present in the k-mer
database but unreachable by the standard walking algorithm.

The split-read k-mer analysis approach:
1. Extract junction k-mers directly from the target definition (not from walking).
2. Query these k-mers against the database.
3. If junction k-mers are present, infer the INDEL.

### Implementation

```
For a large deletion target:
  Reference: ...FLANK_A[DELETED_REGION]FLANK_B...
  Expected mutant: ...FLANK_A|FLANK_B...

  Junction k-mers: All k-mers spanning the junction point.

  Query each junction k-mer against the database.
  If count > threshold: deletion detected.
  Estimate VAF from junction k-mer count / flanking k-mer count.
```

This bypasses the walking algorithm entirely for large INDELs, using direct k-mer
lookup instead of graph-based discovery.

### Advantages

- No k-1 barrier: works for any INDEL size as long as reads span the breakpoint.
- Computationally simple: just k-mer lookups, no graph construction.
- Low false positive rate: junction k-mers are highly specific (novel sequences).

### Limitations

- Requires pre-specified targets: The junction k-mers must be computed from known
  variant definitions. Cannot discover novel large INDELs.
- Limited to INDEL sizes smaller than read_length - k: reads must span the junction.
- Does not produce a graph or path, so detailed variant characterization is limited.
- Cannot distinguish between similar-sized INDELs at nearby positions.

### When to Use

Split-read k-mer analysis is most useful for:
- Known recurrent large INDELs (FLT3-ITD alleles, CALR exon 9 variants).
- Large deletions in tumor suppressor genes (TP53, BRCA1/2).
- Monitoring specific variants identified in a patient's tumor.

## Strategy 3: Tumor-Informed Junction K-mer Search

### Concept

When the patient's tumor genotype is known (from a prior tissue biopsy or initial
diagnostic sample), the exact junction k-mers for each variant can be pre-computed.
Monitoring then reduces to searching for these specific k-mers in each longitudinal
ctDNA sample.

### Workflow

```
1. Tumor profiling (one-time):
   - WGS/WES/panel sequencing of tumor
   - Call variants including large INDELs and SVs
   - For each variant, compute junction k-mers

2. Target generation:
   - Create a k-mer target set: {junction k-mer: variant annotation}
   - Include positive and negative control k-mers

3. Longitudinal monitoring:
   - For each ctDNA sample, create jellyfish database
   - Query all junction k-mers against the database
   - Report: variant, junction k-mer count, estimated VAF

4. MRD assessment:
   - Aggregate across multiple variants for robust MRD call
   - Use statistical model (e.g., Monte Carlo integration) for composite MRD score
```

### Advantages

- Maximum sensitivity: Searches for the exact k-mer signatures of known variants.
- No graph-based discovery needed: Direct lookup is faster and simpler.
- Works for any variant size: Large INDELs, fusions, SVs --- all reduce to junction
  k-mer queries.
- Enables multi-variant MRD scoring: Aggregate signal across 10-50 patient-specific
  variants for robust MRD detection.

### Implementation in kmerdet

This approach maps naturally to a new subcommand:

```bash
kmerdet monitor --junctions patient_junctions.tsv --db ctdna_sample.jf -o mrd_report.tsv
```

Where `patient_junctions.tsv` contains:
```tsv
variant_id          junction_kmer                    ref_kmer                         annotation
FLT3_ITD_32bp       ACGTACGTACGT...TGCATGCA         ACGTACGT...ACGTACGT              FLT3 ITD 32bp
EGFR_del19_15bp     TGCATGCATGCA...ACGTACGT         TGCATGCA...TGCATGCA              EGFR exon 19 del
```

## Strategy 4: Coverage Drop Detection for Large Deletions

### Concept

Large deletions cause a measurable drop in k-mer coverage across the deleted region.
Even when the graph is disconnected (deletion > k-1), the coverage signal is
detectable from existing reference k-mer counts.

### Detection Algorithm

```
For target region with reference k-mers K1, K2, ..., Kn:

1. Compute median coverage: C_median = median(count(Ki) for all i)
2. Compute per-k-mer ratios: Ri = count(Ki) / C_median
3. Identify coverage drops:
   - Sliding window of size W (e.g., 50bp)
   - For each window, compute mean ratio R_window
   - Flag windows where R_window < threshold (e.g., 0.5 for heterozygous deletion)

4. Estimate deletion boundaries:
   - Left boundary: first position where R drops below threshold
   - Right boundary: last position where R drops below threshold
   - Deletion size estimate: right - left

5. Estimate VAF:
   - For heterozygous deletion: VAF ~ 2 * (1 - R_window)
   - For homozygous deletion: VAF ~ 1 - R_window
   - Adjust for tumor fraction
```

### Coverage Drop vs. Normal Variation

The challenge is distinguishing true deletions from normal coverage variation:

```
Normal coverage variation (no deletion):
  Position: 1    2    3    4    5    6    7    8    9    10
  Coverage: 100  95   105  98   102  97   103  99   101  100
  Ratio:    1.0  0.95 1.05 0.98 1.02 0.97 1.03 0.99 1.01 1.0

Heterozygous deletion at 50% tumor fraction (VAF = 25%):
  Position: 1    2    3    4    5    6    7    8    9    10
  Coverage: 100  98   75   74   76   73   77   75   99   101
  Ratio:    1.0  0.98 0.75 0.74 0.76 0.73 0.77 0.75 0.99 1.01
```

At low tumor fraction, the coverage drop is smaller and harder to detect:

| Tumor fraction | Zygosity      | Expected coverage drop | Detectable? |
|----------------|---------------|----------------------|-------------|
| 50%            | Heterozygous  | 25%                  | Yes         |
| 10%            | Heterozygous  | 5%                   | Marginal    |
| 1%             | Heterozygous  | 0.5%                 | No          |
| 50%            | Homozygous    | 50%                  | Yes         |
| 10%            | Homozygous    | 10%                  | Yes         |
| 1%             | Homozygous    | 1%                   | Marginal    |

### Existing Implementation

The `kmerdet coverage` subcommand already reports per-target k-mer coverage statistics.
Extending this to detect coverage drops requires:

1. Per-position coverage output (not just aggregate statistics).
2. Sliding window analysis with configurable window size and threshold.
3. Integration with the main `detect` workflow as a supplementary signal.

### Limitations

- Low sensitivity at low VAF: Coverage drops are proportional to tumor fraction *
  VAF, making them undetectable at the VAF levels relevant to MRD monitoring.
- No breakpoint resolution: Coverage drops indicate the approximate region of a
  deletion but cannot pinpoint exact breakpoints.
- Confounded by GC bias, mappability, and library artifacts.
- Best suited for high-tumor-fraction diagnostic samples, not longitudinal MRD
  monitoring.

## Strategy Integration

The four strategies are complementary and should be applied based on the clinical
context:

```
                          Tumor fraction
                    High (>10%)     Low (<1%)
                  +---------------+---------------+
   Known variant  | Strategy 4    | Strategy 3    |
   (monitoring)   | (cov. drop)   | (junction     |
                  | + Strategy 3  |  k-mer search)|
                  +---------------+---------------+
   Unknown variant| Strategy 1    | Strategy 1    |
   (discovery)    | (multi-k)     | (multi-k)     |
                  | + Strategy 2  | + Strategy 2  |
                  | + Strategy 4  |               |
                  +---------------+---------------+
```

### Recommended Implementation Order

1. **Multi-k approach** (Strategy 1): Highest impact, broadest applicability.
   Already partially implemented. Complete the merge logic and validate with
   benchmarks.

2. **Split-read k-mer analysis** (Strategy 2): Moderate implementation effort,
   high value for known recurrent large INDELs. Natural extension of the existing
   target-based approach.

3. **Tumor-informed junction search** (Strategy 3): Most sensitive for MRD monitoring
   when tumor genotype is known. Requires a new subcommand but conceptually simple.

4. **Coverage drop detection** (Strategy 4): Limited to high-tumor-fraction samples.
   Lowest priority for the liquid biopsy MRD monitoring use case, but useful for
   diagnostic samples.

## Benchmarking Expectations

### INDEL Size vs. Sensitivity Curves

The primary benchmark metric is sensitivity as a function of INDEL size, stratified
by k value and strategy:

```
Target performance (multi-k, k=21 + k=31):

INDEL size (bp):  1-3   4-7   8-15  16-20  21-30  31-50
Sensitivity:      60%   50%   40%   35%    20%    10%

vs. current (k=31 only):
INDEL size (bp):  1-3   4-7   8-15  16-20  21-30  31-50
Sensitivity:      55%   38%   10%   5%     3%     0%
```

### Key Crossover Points

- **k=21 vs. k=31 for INDELs 15-25bp:** k=21 should significantly outperform k=31
  for this size range, as many of these INDELs exceed the k=31 barrier but not the
  k=21 barrier.
- **Multi-k vs. single-k for INDELs >10bp:** Multi-k should show consistent
  improvement across all large INDEL sizes.
- **Specificity crossover:** For INDELs <5bp, k=31 may have higher specificity than
  k=21 due to fewer false positive branches. The merge strategy should preserve this
  advantage.

### Target: >50% Sensitivity for INDELs up to 20bp

With the multi-k approach (k=21 + k=31), the target is to achieve >50% sensitivity
for INDELs up to 20bp at 0.1% VAF with 1000x deduplicated coverage. This would
represent a transformative improvement over the current 3% sensitivity for >7bp
INDELs.

## References

- Baudry, K., et al. (thesis). K-mer-based variant detection for liquid biopsy
  ctDNA monitoring. (Source of sensitivity figures cited in this document.)
- Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences.
  Bioinformatics, 34(18), 3094-3100.
- Rizk, G., Lavenier, D., & Chikhi, R. (2013). DSK: k-mer counting with very low
  memory usage. Bioinformatics, 29(5), 652-653.
- Iqbal, Z., et al. (2012). De novo assembly and genotyping of variants using colored
  de Bruijn graphs. Nature Genetics, 44(2), 226-232.
- Ye, K., et al. (2009). Pindel: a pattern growth approach to detect break points
  of large deletions and medium sized insertions from paired-end short reads.
  Bioinformatics, 25(21), 2865-2871.
- Garrison, E., & Marth, G. (2012). Haplotype-based variant detection from short-read
  sequencing. arXiv preprint arXiv:1207.3907.
