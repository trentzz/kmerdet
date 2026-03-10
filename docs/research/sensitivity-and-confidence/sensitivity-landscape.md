# Sensitivity Landscape: Where Variants Are Being Missed and Why

## Overview

The thesis validation of the kam pipeline on 10 patients with UMI-based duplex sequencing established a baseline sensitivity of 77% for SNVs, substantially lower for short indels (~38%), and near-zero for long indels (~3%). This document dissects the sources of false negatives across variant types, VAF ranges, and pipeline stages, drawing on both the thesis findings and broader literature on k-mer-based and low-VAF variant detection.

---

## 1. Sensitivity Breakdown by Variant Type

### 1.1 SNVs: 77% Sensitivity -- Where Are the Remaining 23% Lost?

For single nucleotide variants, exactly `k` k-mers differ between the reference and variant paths. This is the simplest case for the walking algorithm, yet 23% of true SNVs are still missed. The losses decompose into several distinct failure modes:

#### Very Low VAF (<0.1%)

At VAF 0.1% with typical duplex sequencing depth of 2000-5000x, variant k-mers have expected counts of only 2-5. The km parameter `--count 2` sets the absolute floor, meaning variants producing fewer than 2 supporting k-mers are invisible. At VAF 0.05% with 2000x coverage, the expected variant k-mer count is ~1, which is below the minimum threshold. The binomial probability of observing zero variant k-mers at a true 0.1% VAF with 2000x coverage is non-trivial:

```
P(count = 0) = (1 - 0.001)^2000 ~ 0.135
P(count = 1) = 2000 * 0.001 * (1 - 0.001)^1999 ~ 0.271
P(count < 2) ~ 0.406
```

This means approximately 40% of true variants at 0.1% VAF with 2000x coverage will have k-mer counts below 2 -- a substantial source of false negatives. This aligns with the thesis finding that the practical detection limit is VAF > 0.1%, and even at that boundary, sensitivity drops sharply.

#### High Background Noise Regions

Certain genomic regions produce elevated background k-mer noise due to:

- **High GC content**: Regions with extreme GC composition show uneven coverage, creating pockets where variant k-mers are drowned in noise. GC bias in PCR amplification can cause 10-fold coverage variation across the target.
- **Repetitive elements**: K-mers overlapping SINE/LINE elements, satellite repeats, or segmental duplications may have elevated counts from off-target sources, raising the effective noise floor above the variant signal.
- **Pseudogene interference**: Processed pseudogenes share near-identical sequence with their parent genes. K-mers from pseudogene-containing reads can create spurious counts at the target locus, masking the true variant signal or raising the threshold above the variant's count.

#### Sequencing Error Hotspots Overlapping True Variants

Illumina sequencing exhibits systematic, context-dependent error patterns. The two primary motifs are GGC sequences and inverted repeats (Schirmer et al., NAR 2015; Manley et al., NAR Genomics 2021). When a true variant happens to occur at a position that is also a sequencing error hotspot, the error correction and filtering logic may aggressively remove the variant k-mers along with the error k-mers. Specifically:

- The GGC motif preceding a base position biases the following base toward the same identity, creating systematic C>T and G>A errors at rates 5-10x higher than the background error rate.
- These position-specific error rates can reach 0.5-1%, directly overlapping with the VAF range of interest for liquid biopsy (0.1-1%).
- When error k-mers at a position consistently match a true variant allele, the deduplication and counting pipeline cannot distinguish them.

#### Fixed Threshold Too Aggressive for Borderline Cases

The km walking algorithm uses a relative threshold of `max(sum_of_sibling_counts * ratio, absolute_minimum)` to decide whether to extend a k-mer branch. The fixed `ratio = 0.30` (30% of total sibling count) is calibrated for typical coverage, but fails at the extremes:

- **Low-coverage targets**: When total coverage drops below ~500x, a 0.5% VAF variant has an expected variant k-mer count of ~2.5. The reference k-mer count is ~500. The ratio 2.5/500 = 0.005, far below the 0.30 threshold. The walking algorithm never extends this branch.
- **High-coverage, low-VAF**: At 5000x coverage with 0.1% VAF, the variant count is ~5 while the reference count is ~5000. The ratio 5/5000 = 0.001, again far below threshold.

The key insight: **the walking ratio threshold operates on sibling count ratios, not absolute counts, but for low-VAF detection the ratio will always be far below the 30% default**. The separate `--count` parameter (absolute minimum) is the actual determinant for low-VAF sensitivity, yet it provides no statistical framework for assessing whether the count is significant given the local coverage.

### 1.2 Short Indels (1-7 bp): ~38% Sensitivity

Short indel sensitivity is roughly half that of SNVs. The causes compound:

#### Homopolymer Context Ambiguity

Approximately 20-25% of short indels in the human genome occur in homopolymer runs (Fang et al., 2014). For a 1-bp deletion in a poly-A run of length 10, the resulting variant creates k-mers that may be identical to k-mers from any position within the run. This creates multiple valid representations of the same variant (the "left-alignment ambiguity" problem), and for k-mer based detection, these ambiguous k-mers reduce the effective signal because variant-supporting k-mers cannot be unambiguously assigned to the variant path.

The de Bruijn graph representation further compounds this: in a homopolymer of length L with k-mer length k, if L > k, all k-mers within the run are identical (AAAA...A). A 1-bp deletion changes L to L-1 but does not change any interior k-mers -- only the junction k-mers at the boundaries differ. This means the variant signal comes from just 2 junction k-mers rather than the k junction k-mers expected for an indel in non-repetitive sequence.

Scalpel (Narzisi et al., 2014) demonstrated that assembly-based approaches in micro-assembly windows achieve substantially better indel calling in homopolymer regions compared to alignment-based approaches. The k-mer walking approach, while related to assembly, lacks the localized re-assembly that tools like Scalpel employ.

#### K-mer Walking Branching Explosion

Short indels produce branching points in the k-mer walk where both the reference path and the variant path are supported. For a 3-bp insertion, the walking algorithm encounters:

1. A branch point at the insertion boundary where 4 possible extensions exist
2. Within the insertion, up to 3 additional branch points
3. At the re-convergence point, another branch decision

The total number of paths explored grows combinatorially with the number of branch points. The `max_break = 10` parameter limits branching, but complex indels near other variants or in repetitive contexts can exhaust this budget before the variant path is fully traversed. When the walking algorithm hits `max_stack = 500` or `max_break = 10`, it terminates early, potentially before discovering the variant path.

At low VAF, the problem is worse: the variant branch has much lower k-mer counts than the reference branch, making it more likely to be pruned by the ratio threshold before the walk reaches the re-convergence point.

#### Classification Edge Cases for Complex Indels

The `diff_path_without_overlap` classification logic compares reference and alternative paths to determine variant type. Edge cases include:

- **Tandem repeat expansions/contractions**: A 3-bp deletion in a CAGCAGCAG repeat may be classified as a complex indel rather than a simple deletion because the alignment between paths is ambiguous.
- **Indels near target boundaries**: When the indel occurs within k bases of the target start or end, the variant path may not fully connect to the source/sink anchors, causing classification failure.
- **Multi-nucleotide substitutions adjacent to indels**: An SNV within 5-10 bases of a short indel creates a compound variant that may produce a single path with mixed classification signals.

### 1.3 Long Indels (>7 bp): ~3% Sensitivity

Long indels are nearly undetectable with the current approach, and the causes are fundamental to k-mer-based detection:

#### K-mer Length Constraint

For k=31 and an insertion of length L, the number of k-mers that span the insertion junction is max(0, k - L). For:

| Insertion Length | Spanning K-mers | % of Full Signal (vs. SNV's k=31) |
|-----------------|-----------------|-----------------------------------|
| 5 bp | 26 | 84% |
| 10 bp | 21 | 68% |
| 15 bp | 16 | 52% |
| 20 bp | 11 | 35% |
| 25 bp | 6 | 19% |
| 30 bp | 1 | 3% |
| >30 bp | 0 | 0% -- fully novel k-mers only |

For insertions larger than k (31 bp), zero k-mers span the junction; instead, entirely novel k-mers must be present within the insertion itself. These novel k-mers have no reference counterpart, exist only from the variant allele, and at low VAF produce counts that are easily below threshold.

For deletions, the constraint is different but equally severe: a large deletion removes k-mers from the reference path and creates junction k-mers that span the deletion breakpoint. At low VAF, these junction k-mers may have counts below the walking threshold.

#### Walking Depth Limits

The iterative DFS walking algorithm has a maximum stack depth (`max_stack = 500`). For a long insertion, the walk must traverse:
- Reference k-mers up to the insertion point
- The entire insertion sequence (L additional k-mers)
- Reference k-mers from the insertion endpoint back to the sink

For an insertion of 100 bp, this adds ~100 extra nodes to the walk. Combined with branching at the insertion boundaries and within the insertion itself (where novel k-mers may have ambiguous extensions), the stack depth limit can be reached before the walk completes.

#### Graph Complexity

Long indels produce graphs with many nodes that have limited k-mer count support. The shortest-path algorithm (Dijkstra) may not find the variant path if:
- Any k-mer along the variant path has zero count (path is disconnected)
- The path weight exceeds the maximum search depth
- Alternative shorter paths through the graph are preferred by the cost function

The 2-kupl approach (Pelletier et al., 2021) demonstrated that k-mer-based methods can detect deletions >100 bp by comparing k-mer populations between matched samples rather than walking through the k-mer graph. This suggests that for long indels, a differential k-mer approach may be more effective than the walking approach.

---

## 2. Sensitivity Breakdown by VAF Range

### 2.1 VAF > 1%: Reliable Detection (>90%)

At VAF > 1% with typical duplex sequencing depths (2000-5000x), variant k-mers have counts of 20-50+. The walking algorithm easily extends through variant branches, and the signal-to-noise ratio is favorable. Detection failure at this VAF range is primarily due to:

- Indel-specific issues (as described above) -- the variant type matters more than VAF at this level
- Target design problems: anchor k-mers that are not unique in the genome
- Extreme GC bias causing localized coverage drops
- Rare cases of allelic dropout in library preparation

Benchmarking studies of commercial ctDNA assays show that at VAF > 0.5%, most methods achieve >95% sensitivity for SNVs (Northstar Select validation: 95% LOD at 0.15% VAF for SNVs; Guardant360 and FoundationOne Liquid CDx both report >99% sensitivity above 0.5% VAF).

### 2.2 VAF 0.1-1%: The Critical Detection Zone

This is the clinically most important range for MRD monitoring and early relapse detection. It is also where most false negatives occur in the kam pipeline. The challenges compound:

**Statistical detection limit**: At 0.3% VAF with 3000x coverage, the expected variant k-mer count is ~9. With Poisson-distributed sampling variation, the 95% confidence interval is approximately 3-18. There is a ~3% chance of observing fewer than 3 counts, which would be right at the `--count 2` threshold. At 0.1% VAF, this probability rises to ~40%.

**Signal-to-noise overlap**: The per-base sequencing error rate after UMI deduplication is approximately 0.01-0.1% (depending on duplex vs. simplex consensus). This means error k-mers can have counts in the same range as true variant k-mers at 0.1-0.3% VAF. The pipeline cannot distinguish a true variant k-mer with count 3 from an error k-mer with count 3 based on count alone.

**Coverage non-uniformity**: Even with targeted sequencing, coverage varies across the target region. If the variant falls in a low-coverage pocket within the target, the effective VAF threshold at that position is higher than the nominal 0.1%.

**Modern approaches to improve detection in this range** include:

- **UMI-aware error correction**: Consensus calling within UMI families reduces the effective error rate to <0.001% for duplex consensus, pushing the noise floor well below 0.1% VAF (Schmitt et al., PNAS 2012). The kam pipeline uses HUMID for deduplication but does not perform error correction at the consensus level before k-mer counting.
- **Duplex sequencing**: Requiring support from both strands of a DNA molecule provides error rates of ~10^-7, enabling detection at VAF 0.01% or lower (Kennedy et al., Nature Protocols 2014).
- **CRISPR-based enrichment**: MUTE-Seq uses FnCas9 to selectively remove wild-type DNA, boosting the effective VAF of mutant alleles by 10-100x (demonstrated detection at 0.005% VAF).

### 2.3 VAF < 0.1%: Below Current Detection Limit

At this range, the expected variant k-mer count at typical sequencing depths is 0-2, making reliable detection impossible with standard k-mer counting. What would it take to detect at this level?

**Sequencing depth**: At VAF 0.01%, achieving an expected count of 5 variant k-mers requires 50,000x raw coverage. With typical 150 bp paired-end reads and 167 bp cfDNA fragments, this requires ~7.5 million reads per target region. This is achievable with hybrid-capture targeted sequencing.

**Error suppression**: The background error rate must be pushed below the target VAF. Duplex UMI consensus achieves error rates of ~10^-7, which is sufficient. However, the kam pipeline counts k-mers from deduplicated reads (simplex consensus), not duplex consensus, achieving error rates of only ~10^-4.

**Molecular counting**: Rather than counting k-mer occurrences in reads, counting distinct UMI families supporting each k-mer would provide a molecule-level count that is immune to PCR amplification bias. At 50,000x sequencing depth with 5000 unique molecules, a 0.01% VAF variant would be supported by ~0.5 unique molecules on average -- still borderline.

**Multi-target integration**: By combining evidence across multiple variant targets per patient (e.g., 50 targets), the aggregate signal can provide detection sensitivity below the per-target limit. If each target independently has a 1% chance of detecting the variant, the probability of detecting at least one across 50 targets is 1 - (0.99)^50 ~ 39%. Methods like CAPP-Seq and Phased Variant Enrichment exploit this principle.

---

## 3. Pipeline Stage Analysis: Where Signal Is Lost

### 3.1 Deduplication (HUMID)

**Over-aggressive UMI collapsing**: HUMID identifies UMI families by clustering reads with matching UMI barcodes and similar mapping positions. Potential issues:

- **UMI collisions**: Different original molecules that happen to share the same UMI barcode are collapsed into a single consensus, losing one of the original signals. With 16-bp UMIs (4^16 = ~4 billion possible UMIs), the collision rate is low for typical library sizes. But with 12-bp UMIs (16 million possible) and 5 million unique molecules, the expected collision rate is ~1 per 6.4 unique molecules (by birthday problem approximation), which can cause occasional loss of variant-supporting molecules.
- **Consensus calling errors**: When a UMI family contains a mix of variant-supporting and reference-supporting reads (heterogeneous family), the consensus may call the reference base if the variant reads are in the minority. This systematically reduces the variant signal.
- **Family size threshold**: Requiring a minimum family size for consensus (e.g., >= 3 reads) discards small families that may contain variant-supporting reads. At very low input DNA amounts (common in cfDNA), many families have size 1-2.

**Estimated signal loss at this stage**: 5-15% of variant-supporting molecules may be lost through deduplication, depending on library complexity, UMI design, and family size distribution.

### 3.2 K-mer Counting (Jellyfish)

**Low-count cutoff (`-L 2`)**: Jellyfish's `-L` parameter discards k-mers with counts below the specified threshold. With `-L 2`, any k-mer appearing exactly once is removed from the database.

At very low VAF, a variant k-mer may appear in only 1 read (after deduplication). Setting `-L 2` directly eliminates these singleton k-mers. However, setting `-L 1` would retain enormous numbers of error k-mers (each sequencing error creates a novel k-mer that appears once), dramatically increasing database size and downstream noise.

**The tradeoff**: `-L 2` loses ~30-40% of true variant k-mers at 0.1% VAF (those that happen to have count 1 after deduplication), while `-L 1` would increase false positive k-mers by 10-100x.

**Canonical counting**: Jellyfish's `-C` flag counts k-mers and their reverse complements together. This doubles the effective count for k-mers that are not palindromic (most k-mers), improving sensitivity. However, for palindromic or near-palindromic k-mers, the canonical count may mix signal from different genomic contexts.

### 3.3 Walking (km find_mutation)

**Threshold too high**: The walking extension threshold `max(sibling_sum * ratio, n_cutoff)` determines whether a branch is explored. With `ratio = 0.30`:

- At a branch point with 1000 reference k-mers and 5 variant k-mers, the ratio is 0.005. The ratio test fails (0.005 < 0.30). The absolute minimum `n_cutoff = 2` saves this branch if count >= 2.
- But `n_cutoff` is not the `--count` parameter used in the kam pipeline. In the km source code, `n_cutoff` is derived from `count * ratio`. With `count = 2, ratio = 0.00001`, the effective `n_cutoff = 0.00002`, which is floored to effectively 1 or 2 depending on implementation.

The walking threshold is the single largest source of false negatives for low-VAF variants. When the walk encounters a branch point where the variant extension has a count of 2 and the reference extension has a count of 2000, the variant branch may or may not be explored depending on which threshold dominates (ratio vs. absolute minimum).

**Walking stops before reaching variant**: The DFS walk has bounded depth (`max_stack`). If the target region is long (>300 bp) and the variant is near the center, the walk must traverse ~150 k-mers from each anchor to reach the variant site. With branching from noise, the effective search tree can exhaust the stack limit before reaching the variant position.

### 3.4 Graph Construction

**Variant path not connected source-to-sink**: For the Dijkstra shortest-path algorithm to find the variant path, every k-mer along that path must be present in the graph (count >= threshold). A single missing k-mer creates a gap that disconnects the path.

At low VAF, the probability of any single variant k-mer having count = 0 (after deduplication and counting) is non-negligible. For a variant path of length L k-mers, each with independent probability p of being present:

```
P(path connected) = p^L
```

For L = 31 (typical for an SNV) and p = 0.95: P(connected) = 0.95^31 = 0.21. This means even with 95% per-k-mer detection probability, only 21% of SNV paths are fully connected. In practice, p is higher for SNVs at detectable VAF (p > 0.99), but for indels with more junction k-mers and lower per-k-mer counts, path connectivity is a real concern.

**Edge weight bias**: Reference edges are weighted 0.01 while non-reference edges are weighted 1.0 (100x higher). This strongly biases Dijkstra toward reference paths, which is the intended behavior for finding the shortest reference path. However, it means that variant paths with even one reference-weight edge will be strongly preferred over purely variant paths, potentially causing the algorithm to return chimeric paths that mix reference and variant segments.

### 3.5 Classification (diff_path_without_overlap)

**Edge cases in classification**: The classification logic walks from both ends of the reference and variant paths to find the divergence region. Failures occur when:

- **The variant is at the very start or end of the target**: The "first mismatch from left" search immediately finds the divergence, but there is insufficient flanking sequence to anchor the classification.
- **Multiple variants on the same path**: If two SNVs are within k bases of each other, they produce a single alternative path with two mismatches. The classification logic may report this as a complex indel or MNP rather than two separate SNVs.
- **Variant path has different length than reference but contains substitutions**: A compound indel+SNV event may be misclassified as a simple indel or vice versa, depending on the overlap calculation.

### 3.6 Filtering (kmtools filter)

**True variants filtered as artifacts**: The kmtools filter operates in two modes:

- **Reference mode**: Compares detected variants against a reference database. Variants matching known germline variants or common artifacts are filtered out. True somatic variants that coincidentally match a common polymorphism (e.g., in dbSNP) can be incorrectly filtered.
- **Alt-sequence mode**: Filters based on alternative sequence properties. Variants with rVAF below a threshold, insufficient minimum coverage, or expression below a cutoff are removed.

The thesis used relatively conservative filtering parameters optimized for specificity (zero false positives). This aggressiveness inevitably sacrifices some true positives, particularly:
- True variants at VAF just above the detection threshold
- Variants in regions with higher-than-average background noise
- Variants with borderline expression or coverage metrics

---

## 4. Sensitivity Comparison: K-mer vs. Alignment-Based Approaches

### 4.1 Where K-mer Methods Win

- **Complex structural variants**: ITDs, gene fusions, and complex rearrangements that produce chimeric reads are often handled better by k-mer walking, which naturally discovers novel paths without requiring proper alignment of chimeric fragments.
- **Repetitive regions**: K-mers spanning a variant in a repetitive region are retained regardless of mapping quality. Alignment-based callers may filter these as low-MAPQ.
- **Speed**: The pre-computed k-mer database enables targeted queries in seconds. km processes a 50-target panel in ~2 minutes vs. 30-60 minutes for GATK HaplotypeCaller.

### 4.2 Where Alignment-Based Methods Win

- **Indels**: Assembly-based callers like Scalpel, Lancet, and GATK HaplotypeCaller (which performs local reassembly) consistently outperform k-mer walking for indels, especially in repetitive contexts. Benchmarks show alignment-based indel sensitivity of 80-95% vs. ~38% for k-mer walking.
- **Low VAF SNVs**: Specialized somatic callers like Mutect2, Strelka2, and VarDict use sophisticated statistical models (Bayesian, Fisher's exact test, strand bias) that are absent from the k-mer counting approach. At VAF 0.1-0.5%, these callers achieve 80-95% sensitivity vs. ~70% for the k-mer approach.
- **Quality metrics**: Alignment-based callers produce QUAL, GQ, MAPQ, strand bias, and other quality metrics that enable post-hoc filtering. The k-mer approach provides only rVAF, expression, and min_coverage, limiting the ability to distinguish borderline true positives from artifacts.

### 4.3 Comparative Sensitivity Data

| Variant Type | K-mer (km/kam) | Mutect2 | Strelka2 | VarDict | GATK HC |
|-------------|---------------|---------|----------|---------|---------|
| SNVs (>1% VAF) | ~95% | 98% | 97% | 96% | 95% |
| SNVs (0.1-1% VAF) | ~60% | 85% | 80% | 88% | 70% |
| Short indels (1-7 bp) | ~38% | 80% | 78% | 75% | 82% |
| Long indels (>7 bp) | ~3% | 55% | 50% | 60% | 65% |

*Values for k-mer approach from thesis validation; alignment-based values estimated from published benchmarks (Barbitoff et al., 2022; Nature Communications 2025 ctDNA benchmark). Exact comparisons are difficult due to different cohorts and target panels.*

Recent comprehensive benchmarking by Nature Communications (2025) across colorectal and breast cancer liquid biopsy samples confirmed that mutations at VAF < 0.5% remain challenging for all methods, with sensitivity varying significantly across tools and variant types.

---

## 5. Strategies for Improving Sensitivity

### 5.1 Immediate Improvements (kmerdet Phase 1-3)

1. **Adaptive thresholds**: Replace fixed `ratio` and `count` with coverage-dependent thresholds. At each branch point, compute a binomial p-value for the observed variant k-mer count given the local coverage and expected error rate. Extend branches with p < 0.01 regardless of ratio.

2. **Per-target k optimization**: Use shorter k (21-25) for targets with known large indels. The thesis recommends this explicitly (Chapter 6.5). This directly addresses the k-mer length constraint for long indels.

3. **Multi-pass detection**: Run the walking algorithm with multiple parameter sets (permissive first, then stringent) and take the union of detected variants with appropriate false discovery rate control.

### 5.2 Medium-Term Improvements (kmerdet Phase 4-6)

4. **Duplex consensus k-mer counting**: Count k-mers from duplex consensus sequences rather than deduplicated reads, reducing the effective error rate from ~10^-4 to ~10^-7 and pushing the detection limit to VAF ~0.01%.

5. **Junction k-mer enrichment for indels**: Instead of walking through the entire indel, specifically search for junction k-mers that span the insertion/deletion breakpoint. This requires knowing the expected variant sequence (tumor-informed approach) but is far more sensitive than discovery-based walking.

6. **Statistical walking threshold**: Replace the fixed ratio threshold with a likelihood ratio test comparing "variant branch" vs. "noise branch" hypotheses, incorporating the local coverage, error rate, and VAF prior.

### 5.3 Long-Term Improvements (kmerdet Phase 7+)

7. **UMI-aware molecular counting**: Count distinct molecules (UMI families) supporting each k-mer rather than total reads. This eliminates PCR amplification bias and provides a more accurate estimate of the true variant molecule count.

8. **Multi-target integration**: Combine evidence across all targets for a patient to compute a patient-level detection score. Even if no single target achieves significance, the aggregate signal across 50 targets may be significant.

9. **Machine learning variant scorer**: Train a gradient-boosted classifier on features including k-mer count, rVAF, coverage uniformity, sequence context, GC content, and homopolymer proximity to distinguish true variants from artifacts without hard thresholds.

---

## References

- Audoux et al. "DE-kupl: exhaustive capture of biological variation in RNA-seq data through k-mer decomposition." Genome Biology, 2017.
- Pelletier et al. "[2-kupl: mapping-free variant detection from DNA-seq data of matched samples](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04185-6)." BMC Bioinformatics, 2021.
- Barbitoff et al. "[Systematic benchmark of state-of-the-art variant calling pipelines](https://link.springer.com/article/10.1186/s12864-022-08365-3)." BMC Genomics, 2022.
- Deng et al. "[Benchmarking UMI-aware and standard variant callers for low frequency ctDNA variant detection](https://link.springer.com/article/10.1186/s12864-024-10737-w)." BMC Genomics, 2024.
- "[Comprehensive benchmarking of methods for mutation calling in circulating tumor DNA](https://www.nature.com/articles/s41467-025-67842-x)." Nature Communications, 2025.
- Zhang et al. "[Improvement of the sensitivity of circulating tumor DNA-based liquid biopsy](https://www.explorationpub.com/Journals/etat/Article/1002333)." Exploration of Targeted Anti-tumor Therapy, 2024.
- "[GeneBits: ultra-sensitive tumour-informed ctDNA monitoring](https://link.springer.com/article/10.1186/s12967-025-06993-3)." Journal of Translational Medicine, 2025.
- "[Validation of a liquid biopsy assay with increased sensitivity](https://www.sciencedirect.com/science/article/pii/S2950195425000384)." 2025.
- Narzisi et al. "[Indel variant analysis of short-read sequencing data with Scalpel](https://pmc.ncbi.nlm.nih.gov/articles/PMC5507611/)." Current Protocols in Bioinformatics, 2018.
- "[Reducing INDEL calling errors in whole genome and exome sequencing data](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-014-0089-z)." Genome Medicine, 2014.
- "[The Challenge of Small-Scale Repeats for Indel Discovery](https://frontiersin.org/articles/10.3389/fbioe.2015.00008/full)." Frontiers in Bioengineering and Biotechnology, 2015.
- "[Analytical evaluation of circulating tumor DNA sequencing assays](https://www.nature.com/articles/s41598-024-54361-w)." Scientific Reports, 2024.
- Narzisi et al. "Lancet: Microassembly based somatic variant caller." [GitHub](https://github.com/nygenome/lancet).
- Thesis Chapters 4-6: Validation Study Design, Results, and Discussion.
