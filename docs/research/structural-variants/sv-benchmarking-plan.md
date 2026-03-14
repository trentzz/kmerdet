# Benchmarking Plan for Structural Variant Detection

## Overview

This document defines a comprehensive benchmarking plan for evaluating kmerdet's
structural variant detection capabilities. The plan covers simulation strategies for
each SV type, comparison methodologies, ground truth formats, and integration with the
existing benchmarking framework.

All benchmarks should be reproducible, automated, and extensible as new SV detection
strategies are implemented. Results should feed directly into development decisions
about which strategies to prioritize.

## General Benchmarking Principles

### Controlled Variables

Every benchmark must control for:
- **VAF:** The variant allele frequency, simulated at 5%, 1%, 0.5%, 0.1%, 0.05%,
  and 0.01%.
- **Coverage:** Total sequencing depth, tested at 100x, 500x, 1000x, 5000x, and
  10000x (post-deduplication).
- **K-mer length:** k=21, k=31, k=41, and multi-k combinations.
- **Error rate:** Base-call error rate in simulated reads (0.1%, 0.5%, 1.0%).
- **Fragment length:** cfDNA fragment length distribution (mean 170bp, SD 20bp for
  liquid biopsy simulation).

### Metrics

For each benchmark, compute:
- **Sensitivity (recall):** Fraction of true variants detected.
- **Specificity:** Fraction of negative targets correctly reported as negative.
- **Precision (PPV):** Fraction of detected variants that are true positives.
- **F1 score:** Harmonic mean of sensitivity and precision.
- **VAF accuracy:** Correlation and RMSE between true VAF and estimated rVAF.
- **Size accuracy:** For INDELs and ITDs, correlation between true size and estimated
  size.
- **Breakpoint accuracy:** For fusions and large deletions, distance between true and
  estimated breakpoint position.

### Statistical Rigor

- Each condition is replicated 10 times with different random seeds.
- Confidence intervals (95%) are reported for all metrics.
- Power analysis determines minimum replicate count for detecting a 10% sensitivity
  difference.
- Multiple testing correction (Bonferroni or Benjamini-Hochberg) for comparisons
  across conditions.

## Benchmark 1: Gene Fusion Simulation

### Synthetic Chimeric Sequence Generation

Create chimeric sequences from known gene fusion pairs. Use real exon sequences from
the GRCh38 reference genome.

**Fusion pairs to simulate:**

| Fusion           | Gene A region        | Gene B region        | Breakpoint variants |
|------------------|---------------------|---------------------|---------------------|
| BCR-ABL1         | BCR exon 13, 14, 1  | ABL1 exon 2          | e13a2, e14a2, e1a2  |
| EML4-ALK         | EML4 exon 6, 13, 20 | ALK exon 20          | v1, v2, v3a         |
| PML-RARA         | PML exon 3, 6       | RARA exon 3          | bcr1, bcr3          |
| ROS1-CD74        | CD74 exon 6         | ROS1 exon 34         | Standard            |
| KIF5B-RET        | KIF5B exon 15, 24   | RET exon 12          | Standard            |

**Simulation procedure:**

```python
def simulate_fusion(gene_a_seq, gene_b_seq, breakpoint_a, breakpoint_b,
                    vaf, coverage, read_length=150, fragment_mean=170):
    """
    Generate a mixed read set containing fusion and normal reads.

    1. Create chimeric sequence: gene_a_seq[:breakpoint_a] + gene_b_seq[breakpoint_b:]
    2. Generate 'coverage * vaf' reads from chimeric sequence
    3. Generate 'coverage * (1 - vaf)' reads from normal gene A and gene B
    4. Add base-call errors at specified rate
    5. Trim to fragment length distribution
    6. Output: FASTA/FASTQ reads + ground truth junction positions
    """
```

**Read tiling across the junction:**

Reads are positioned such that a controlled fraction spans the fusion junction:

```
Fragment positions relative to junction (J):
  [-----read-----]                       (entirely in gene A, no junction info)
           [-----read-----]              (spans junction, contains junction k-mers)
                    [-----read-----]     (spans junction, contains junction k-mers)
                             [-----read-----]  (entirely in gene B, no junction info)
```

At fragment_mean=170bp and breakpoint at the center of a 300bp target, approximately
60-80% of fragments from the fusion allele will span the junction.

**Control targets:**
- Normal (non-fusion) targets from the same genes: gene A reference region, gene B
  reference region.
- These should always report as reference-only (no variant detected).

**Metrics specific to fusion benchmarks:**
- Sensitivity: Fraction of simulated fusions detected at each VAF level.
- Breakpoint position accuracy: Distance (in bp) between true and called breakpoint.
- Partner gene identification accuracy: Are both genes correctly identified?
- Junction k-mer count: Number of junction k-mers detected vs. expected.

### Reciprocal Fusion Testing

For each fusion, also simulate the reciprocal product:

```
Forward:    BCR[1:breakpoint] + ABL1[breakpoint:]
Reciprocal: ABL1[1:breakpoint] + BCR[breakpoint:]
```

Test whether both orientations are independently detectable and whether detecting
both increases confidence.

### Expected Results

| VAF   | Coverage | Expected sensitivity | Notes                              |
|-------|----------|---------------------|------------------------------------|
| 5%    | 1000x   | >95%                | Clear junction k-mer signal        |
| 1%    | 1000x   | >90%                | ~10 junction k-mers expected       |
| 0.1%  | 1000x   | ~50%                | ~1 junction k-mer, at threshold    |
| 0.01% | 1000x   | <10%                | Below detection limit              |
| 0.1%  | 10000x  | >80%                | ~10 junction k-mers at high depth  |
| 0.01% | 10000x  | ~40%                | ~1 junction k-mer                  |

## Benchmark 2: Large Deletion Simulation

### Deletion Size Series

Simulate deletions of increasing size to map the sensitivity curve:

**Deletion sizes:** 5, 10, 15, 20, 25, 30, 35, 40, 50, 75, 100, 150, 200, 500bp

**Target construction:**

```python
def simulate_large_deletion(reference_seq, del_start, del_size,
                             vaf, coverage, read_length=150):
    """
    1. Create mutant sequence: ref[:del_start] + ref[del_start+del_size:]
    2. Generate reads from both reference (1-vaf) and mutant (vaf) sequences
    3. Create jellyfish database from combined reads
    4. Create target FASTA spanning the deletion region
    5. Ground truth: deletion at del_start, size del_size
    """
```

**Reference sequences for deletion simulation:**
- Use 1000bp regions from well-characterized genes (TP53, BRCA1, APC, EGFR).
- Avoid repetitive regions (Alu, LINE) to isolate the effect of deletion size.
- Include some targets in repetitive regions to test robustness.

### K-mer Length Comparison

For each deletion size and VAF, run detection at:
- k=21 only
- k=31 only
- k=41 only
- k=21 + k=31 (multi-k merge)
- k=21 + k=31 + k=41 (triple-k merge)

**Expected sensitivity patterns:**

```
Deletion    k=21   k=31   k=41   k=21+31  k=21+31+41
5bp         80%    85%    85%    88%      88%
10bp        70%    75%    80%    80%      82%
15bp        55%    65%    70%    70%      72%
20bp        40%    55%    65%    58%      66%
25bp        10%    40%    55%    42%      56%
30bp        0%     20%    45%    22%      46%
35bp        0%     5%     30%    7%       32%
50bp        0%     0%     5%     0%       7%
100bp       0%     0%     0%     0%       0%
```

Note: Sensitivity estimates assume walking-based detection. Split-read k-mer analysis
(Strategy 2 from large-indel-strategies.md) would recover some of the 0% cases.

### Graph Connectivity Analysis

For each deletion size and k value, measure:
- **Graph connectivity:** Is the graph connected (single component) or disconnected?
- **Path count:** Number of alternative paths found.
- **Junction k-mer presence:** Are junction k-mers present in the database?
- **Junction k-mer reachability:** Are junction k-mers reachable from reference k-mers
  via walking?

This analysis directly tests the k-1 barrier prediction:

```
k=31: Graph should disconnect at deletion size > 30bp
k=21: Graph should disconnect at deletion size > 20bp
```

### Coverage Drop Detection

For deletions >50bp, additionally evaluate the coverage drop detection strategy:

```
For each deletion:
  1. Compute per-position reference k-mer coverage
  2. Apply sliding window analysis
  3. Report: detected (Y/N), estimated deletion region, estimated size
  4. Compare with true deletion region
```

Metrics:
- Detection rate: Fraction of large deletions flagged by coverage drop.
- Size accuracy: RMSE of estimated vs. true deletion size.
- Boundary accuracy: Distance between estimated and true deletion boundaries.

## Benchmark 3: Large Insertion Simulation

### Insertion Size Series

**Insertion sizes:** 5, 10, 15, 20, 25, 30, 35, 40, 50, 75, 100bp

**Insertion sequences:**
- Random sequence (no homology to reference).
- Partially homologous sequence (50% identity to adjacent reference).
- Repetitive sequence (tandem repeat units).
- Clinically relevant: CALR exon 9 insertions (5bp frameshift), NPM1 exon 12
  insertions (4bp TCTG), FLT3-ITD-like tandem duplications.

**Target construction:**

```python
def simulate_large_insertion(reference_seq, ins_position, ins_sequence,
                              vaf, coverage, read_length=150):
    """
    1. Create mutant: ref[:ins_position] + ins_sequence + ref[ins_position:]
    2. Generate mixed reads at specified VAF
    3. Create jellyfish database
    4. Create target FASTA spanning insertion site
    5. Ground truth: insertion at ins_position, sequence ins_sequence
    """
```

### Metrics Specific to Insertions

- **Detection rate:** Was the insertion detected at all?
- **Size accuracy:** |estimated_size - true_size| in bp.
- **Sequence accuracy:** Edit distance between estimated and true inserted sequence
  (when the algorithm reports the inserted sequence).
- **rVAF accuracy:** |estimated_rVAF - true_VAF|.

### Expected Challenges

```
Insertion size vs. detection challenge:

1-20bp (< k-1 for k=21):
  - Fully bridgeable at k=21
  - Some disconnection at k=31 for 25-30bp
  - Expected sensitivity: 40-80% depending on VAF

21-30bp (< k-1 for k=31):
  - Fully bridgeable at k=31
  - Disconnected at k=21
  - Expected sensitivity: 20-60% at k=31

31-50bp (> k-1 for k=31):
  - Disconnected at all standard k values
  - Requires split-read k-mer analysis
  - Expected sensitivity with walking alone: <5%
  - Expected sensitivity with split-read: 30-50%

>50bp:
  - Disconnected at all standard k values
  - Some insertions exceed read length (no junction reads)
  - Expected sensitivity: <10% even with split-read
```

## Benchmark 4: ITD Simulation

### ITD Size Series

Internal tandem duplications are modeled by duplicating a region of the reference
sequence and inserting it adjacent to the original:

**ITD sizes:** 3, 5, 10, 15, 20, 25, 30, 35, 40, 50, 75, 100bp

**Target gene:** FLT3 exon 14-15 (the canonical ITD region). Use the actual FLT3
reference sequence for biological relevance.

**Simulation:**

```python
def simulate_itd(reference_seq, dup_start, dup_size, vaf, coverage):
    """
    1. Extract duplicated region: dup_seq = ref[dup_start:dup_start+dup_size]
    2. Create mutant: ref[:dup_start+dup_size] + dup_seq + ref[dup_start+dup_size:]
    3. Generate mixed reads
    4. Create jellyfish database
    5. Ground truth: ITD at dup_start, size dup_size
    """
```

### ITD-Specific Metrics

- **Detection rate:** Was the ITD detected?
- **Classification accuracy:** Was it classified as ITD (not insertion or complex)?
- **Size accuracy:** |estimated_size - true_size|.
- **rVAF accuracy:** |estimated_rVAF - true_VAF|.

### Comparison with km ITD Detection

The thesis specifically validates FLT3-ITD detection. Benchmark against km's
performance on the same simulated data:

```bash
# km detection
km find_mutation counts.jf target.fa > km_results.tsv

# kmerdet detection
kmerdet detect --db counts.jf --targets target.fa -o kmerdet_results.tsv

# Compare
diff km_results.tsv kmerdet_results.tsv
```

For each ITD size and VAF, compare:
- Detection rate (km vs. kmerdet)
- Size accuracy (km vs. kmerdet)
- rVAF accuracy (km vs. kmerdet)
- Runtime (km vs. kmerdet)

### ITD Size Heterogeneity

Simulate samples with multiple ITD alleles:

```
Sample 1: Single ITD, 15bp, VAF 5%
Sample 2: Single ITD, 30bp, VAF 1%
Sample 3: Two ITDs, 15bp (VAF 3%) + 30bp (VAF 2%)
Sample 4: Three ITDs, 10bp (VAF 2%) + 20bp (VAF 1.5%) + 40bp (VAF 0.5%)
```

Test whether kmerdet can:
- Detect all ITD alleles in a polyclonal sample.
- Correctly estimate individual allele VAFs.
- Distinguish ITD alleles of different sizes.

## Benchmark 5: Comparison with Dedicated SV Tools

### Tools for Comparison

| Tool    | Version | Input    | SV Types              | Reference                |
|---------|---------|----------|-----------------------|--------------------------|
| Manta   | 1.6+    | BAM      | DEL, INS, DUP, INV, BND | Chen et al. 2016       |
| Delly   | 1.1+    | BAM      | DEL, DUP, INV, TRA, INS | Rausch et al. 2012     |
| GRIDSS  | 2.13+   | BAM      | All SV types          | Cameron et al. 2017      |
| Pindel  | 0.2.5+  | BAM      | Large DEL, INS, INV, TD | Ye et al. 2009          |
| SvABA   | 1.1+    | BAM      | All SV types          | Wala et al. 2018         |

### Comparison Design

These tools operate on BAM files, not k-mer databases, so the comparison requires
generating BAM files from the same simulated reads used for k-mer databases.

**Pipeline:**

```bash
# 1. Simulate reads
python simulate_sv_reads.py --sv-type deletion --size 50 --vaf 0.01 \
    --coverage 1000 --output reads.fq

# 2. Create jellyfish database
jellyfish count -m 31 -s 1G -C reads.fq -o counts.jf

# 3. Align reads (for BAM-based tools)
bwa mem -t 8 reference.fa reads.fq | samtools sort -o aligned.bam
samtools index aligned.bam

# 4. Run kmerdet
kmerdet detect --db counts.jf --targets targets/ -o kmerdet.tsv

# 5. Run comparison tools
manta --bam aligned.bam --reference reference.fa --output manta/
delly call -g reference.fa aligned.bam -o delly.bcf
gridss -r reference.fa -o gridss.vcf aligned.bam

# 6. Compare results
python compare_sv_results.py --truth ground_truth.vcf \
    --kmerdet kmerdet.tsv --manta manta/results.vcf \
    --delly delly.bcf --gridss gridss.vcf
```

### What We Are Comparing

The comparison is NOT intended to show that kmerdet is better than dedicated SV
callers (it likely is not, for large SVs). The comparison is intended to:

1. **Quantify the gap:** How much worse is k-mer walking for each SV type/size
   compared to dedicated tools?
2. **Identify sweet spots:** Where does kmerdet perform comparably (small-medium
   INDELs, fusions with appropriate targets)?
3. **Assess complementarity:** Do kmerdet and dedicated tools detect different
   subsets of variants? Would combining them improve overall sensitivity?
4. **Validate the niche:** Confirm that kmerdet's advantage is in targeted,
   high-sensitivity, low-VAF detection for known variant sites, while dedicated
   tools are better for discovery of novel SVs.

### Expected Outcomes

```
                    Small INDEL  Large INDEL  Fusion     Large DEL
                    (1-20bp)     (20-100bp)   (targeted) (>100bp)
kmerdet (k=31)      Good         Poor         Good       Poor
kmerdet (multi-k)   Good         Fair         Good       Poor
Manta               Fair         Good         N/A        Excellent
Delly               Fair         Good         N/A        Excellent
GRIDSS              Fair         Excellent    N/A        Excellent
Pindel              Good         Good         N/A        Good

Key: Excellent (>90%), Good (60-90%), Fair (30-60%), Poor (<30%)
```

## Ground Truth Formats

### VCF 4.3 for SVs

```vcf
##fileformat=VCFv4.3
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of mate breakend">
##INFO=<ID=GENE_A,Number=1,Type=String,Description="Upstream fusion partner">
##INFO=<ID=GENE_B,Number=1,Type=String,Description="Downstream fusion partner">
##INFO=<ID=TRUE_VAF,Number=1,Type=Float,Description="Simulated true VAF">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE
chr13   28034105  FLT3_ITD_30bp  N  <DUP:TANDEM>  .  .  SVTYPE=DUP;SVLEN=30;END=28034135;TRUE_VAF=0.01  GT  0/1
chr9    5073770   BCR_ABL1_e13a2  N  [chr22:23524426[N  .  .  SVTYPE=BND;MATEID=ABL1_BCR;GENE_A=BCR;GENE_B=ABL1;TRUE_VAF=0.005  GT  0/1
chr7    55174014  EGFR_del19  ATCTCCGAAAGCCAACAAGGAAATCCTCG  A  .  .  SVTYPE=DEL;SVLEN=-28;TRUE_VAF=0.02  GT  0/1
```

### BED Format for Deletion Regions

```bed
#chrom  start       end         name              score  strand  true_vaf  del_size
chr17   7577120     7577220     TP53_del_100bp    0      .       0.01      100
chr13   32936830    32936930    BRCA2_del_100bp   0      .       0.005     100
chr7    55174014    55174042    EGFR_del19_28bp   0      .       0.02      28
```

### Custom TSV Matching kmerdet Output

```tsv
target_name         variant_type  variant_name       ref_expression  alt_expression  rVAF    min_coverage  true_vaf  true_size
FLT3_ITD_target     ITD           FLT3_ITD_30bp      0.99            0.01            0.01    500           0.01      30
BCR_ABL1_e13a2      Fusion        BCR-ABL1_e13a2     0.995           0.005           0.005   480           0.005     N/A
EGFR_del19_target   Deletion      EGFR_p.E746_A750del 0.98           0.02            0.02    520           0.02      28
TP53_ref_control    Reference     TP53_reference      1.0            0.0             0.0     510           0.0       N/A
```

## Integration with Existing Benchmarking Framework

### New Benchmark Preset

Add to `docs/benchmarking/datasets/create_benchmark_data.py`:

```python
class SVBenchmarkPreset:
    """Generate SV benchmark datasets."""

    def __init__(self):
        self.sv_types = ["fusion", "large_deletion", "large_insertion", "itd"]
        self.sizes = {
            "fusion": [None],  # Size not applicable
            "large_deletion": [10, 20, 30, 50, 75, 100, 200, 500],
            "large_insertion": [10, 20, 30, 50, 75, 100],
            "itd": [5, 10, 15, 20, 30, 50, 75, 100],
        }
        self.vafs = [0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001]
        self.coverages = [100, 500, 1000, 5000, 10000]
        self.k_values = [21, 31, 41]

    def generate_all(self, output_dir):
        for sv_type in self.sv_types:
            for size in self.sizes[sv_type]:
                for vaf in self.vafs:
                    for coverage in self.coverages:
                        self.generate_one(output_dir, sv_type, size, vaf, coverage)

    def generate_one(self, output_dir, sv_type, size, vaf, coverage):
        # Generate simulated reads
        # Create jellyfish databases at each k value
        # Create target FASTAs
        # Write ground truth files
        pass
```

### Benchmark Configuration

Add to `docs/benchmarking/framework/benchmark_config.toml`:

```toml
[sv_benchmark]
enabled = true
output_dir = "results/sv_benchmark"

[sv_benchmark.fusion]
gene_pairs = [
    { gene_a = "BCR", gene_b = "ABL1", breakpoints = ["e13a2", "e14a2"] },
    { gene_a = "EML4", gene_b = "ALK", breakpoints = ["v1", "v2"] },
]
vafs = [0.05, 0.01, 0.005, 0.001]
coverages = [1000, 5000, 10000]

[sv_benchmark.large_deletion]
sizes = [10, 20, 30, 50, 100, 200]
vafs = [0.05, 0.01, 0.005, 0.001]
coverages = [1000, 5000]
k_values = [21, 31, 41]

[sv_benchmark.itd]
sizes = [5, 10, 15, 20, 30, 50]
vafs = [0.05, 0.01, 0.005, 0.001]
coverages = [1000, 5000]
k_values = [21, 31, 41]
reference_gene = "FLT3"
```

### Results Analysis Extension

Add SV-specific analysis to the results processing pipeline:

```python
class SVBenchmarkAnalyzer:
    """Analyze SV benchmark results."""

    def sensitivity_vs_size(self, results, sv_type, k_value, vaf):
        """Plot sensitivity as a function of SV size."""
        pass

    def sensitivity_vs_vaf(self, results, sv_type, size, k_value):
        """Plot sensitivity as a function of VAF for a given SV size."""
        pass

    def multi_k_improvement(self, results, sv_type, vaf):
        """Compare single-k vs multi-k sensitivity."""
        pass

    def breakpoint_accuracy(self, results, sv_type):
        """Plot breakpoint position error distribution."""
        pass

    def vaf_accuracy(self, results, sv_type):
        """Plot estimated vs true VAF correlation."""
        pass

    def graph_connectivity(self, results, sv_type, k_value):
        """Plot graph connectivity vs SV size."""
        pass
```

### SV-Specific Plots

The benchmarking framework should produce the following plots:

**Plot 1: Sensitivity vs. SV Size (by k value)**
```
X-axis: SV size (bp), log scale
Y-axis: Sensitivity (0-100%)
Lines: k=21, k=31, k=41, multi-k
Facets: SV type (deletion, insertion, ITD)
Fixed: VAF = 1%, coverage = 1000x
```

**Plot 2: Sensitivity vs. VAF (by SV size)**
```
X-axis: VAF (%), log scale
Y-axis: Sensitivity (0-100%)
Lines: SV size categories (small, medium, large)
Facets: SV type
Fixed: k=31, coverage = 1000x
```

**Plot 3: Multi-k Improvement Heatmap**
```
X-axis: SV size (bp)
Y-axis: VAF (%)
Color: Sensitivity improvement (multi-k vs. best single-k)
Facets: SV type
```

**Plot 4: Fusion Detection Rate vs. Coverage**
```
X-axis: Coverage (x)
Y-axis: Sensitivity (0-100%)
Lines: Fusion pairs (BCR-ABL, EML4-ALK, etc.)
Fixed: VAF = 0.1%
```

**Plot 5: Breakpoint Accuracy Distribution**
```
Histogram of |estimated_breakpoint - true_breakpoint| in bp
Facets: SV type, k value
```

**Plot 6: VAF Estimation Accuracy**
```
Scatter plot: True VAF (x) vs. Estimated rVAF (y)
Identity line for reference
Color: SV type
Facets: k value
```

**Plot 7: Tool Comparison (kmerdet vs. dedicated SV callers)**
```
Grouped bar chart: Sensitivity by SV type and size category
Groups: kmerdet (k=31), kmerdet (multi-k), Manta, Delly, GRIDSS
```

## Execution Plan

### Phase 1: Simulation Infrastructure (Week 1-2)

1. Implement read simulation for each SV type in Python.
2. Validate simulated reads produce expected k-mer patterns.
3. Create jellyfish database generation pipeline.
4. Create target FASTA generation for each SV type.
5. Implement ground truth file generation in all three formats (VCF, BED, TSV).

### Phase 2: Baseline Benchmarks (Week 2-3)

1. Run kmerdet at k=31 on all SV types and sizes.
2. Establish baseline sensitivity curves.
3. Verify that small INDELs (<30bp) match expected sensitivity from the thesis.
4. Identify the exact k-1 barrier cutoff empirically.
5. Document baseline results.

### Phase 3: Multi-k Benchmarks (Week 3-4)

1. Generate jellyfish databases at k=21, k=31, k=41.
2. Run kmerdet at each k value independently.
3. Run multi-k merge.
4. Quantify improvement from multi-k at each SV size.
5. Identify optimal k-value combinations.

### Phase 4: Strategy-Specific Benchmarks (Week 4-6)

1. Implement and benchmark split-read k-mer analysis for large INDELs.
2. Implement and benchmark tumor-informed junction search.
3. Implement and benchmark coverage drop detection.
4. Compare strategy performance across SV types and sizes.

### Phase 5: Tool Comparison (Week 6-8)

1. Set up Manta, Delly, GRIDSS, Pindel pipelines.
2. Run all tools on the same simulated datasets.
3. Compare sensitivity, specificity, and runtime.
4. Identify complementary strengths.
5. Write comparison report.

### Phase 6: Reporting and Documentation (Week 8-9)

1. Generate all plots.
2. Write benchmark report with key findings.
3. Update development roadmap based on benchmark results.
4. Identify highest-impact improvements to prioritize.

## Resource Requirements

### Compute

- Simulation and jellyfish counting: ~2 CPU-hours per condition.
- kmerdet detection: ~0.1 CPU-hours per condition.
- Dedicated SV tools (Manta, GRIDSS, etc.): ~1-4 CPU-hours per condition.
- Total estimate: ~500-1000 CPU-hours for full benchmark suite.
- Recommended: Run on 16-32 core machine, complete in ~1-3 days.

### Storage

- Simulated reads: ~1-10GB per condition (depending on coverage).
- Jellyfish databases: ~100MB-1GB per condition per k value.
- BAM files: ~1-10GB per condition.
- Total estimate: ~500GB-1TB for full benchmark suite.
- Recommendation: Use temporary storage; keep only results and ground truth.

### Software Dependencies

- Python 3.8+ with numpy, scipy, matplotlib, seaborn, pandas.
- Rust toolchain for kmerdet compilation.
- Jellyfish 2.x for k-mer counting.
- BWA + samtools for read alignment (tool comparison).
- Manta, Delly, GRIDSS, Pindel (tool comparison only).
- Snakemake or Nextflow for pipeline orchestration (recommended).

## References

- Chen, X., et al. (2016). Manta: rapid detection of structural variants and indels
  for germline and cancer sequencing applications. Bioinformatics, 32(8), 1220-1222.
- Rausch, T., et al. (2012). DELLY: structural variant discovery by integrated
  paired-end and split-read analysis. Bioinformatics, 28(18), i333-i339.
- Cameron, D.L., et al. (2017). GRIDSS: sensitive and specific genomic rearrangement
  detection using positional de Bruijn graph assembly. Genome Research, 27(12),
  2050-2060.
- Ye, K., et al. (2009). Pindel: a pattern growth approach to detect break points of
  large deletions and medium sized insertions from paired-end short reads.
  Bioinformatics, 25(21), 2865-2871.
- Wala, J.A., et al. (2018). SvABA: genome-wide detection of structural variants and
  indels by local assembly. Genome Research, 28(4), 581-591.
- Zook, J.M., et al. (2020). A robust benchmark for detection of germline large
  deletions and insertions. Nature Biotechnology, 38(11), 1347-1355.
- Chaisson, M.J.P., et al. (2019). Multi-platform discovery of haplotype-resolved
  structural variation in human genomes. Nature Communications, 10(1), 1784.
