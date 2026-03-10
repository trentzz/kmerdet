# K-mer Length Tradeoffs for Variant Detection

## Theoretical Foundation

### K-mer Uniqueness and Specificity

The probability that a specific k-mer occurs by chance in a genome of size G follows from the assumption that the genome is a random sequence over a 4-letter alphabet:

```
P(k-mer occurs at a given position) = 1 / 4^k
P(k-mer appears somewhere in genome) ~= 1 - (1 - 1/4^k)^G ~= G / 4^k   (for G << 4^k)
```

For the human genome (G ~ 3.1 billion):

| k   | 4^k            | P(random match)  | Practical uniqueness |
|-----|----------------|-------------------|----------------------|
| 15  | 1.07 x 10^9    | ~2.9 (expect ~3 matches) | Not unique |
| 17  | 1.72 x 10^10   | ~0.18             | Mostly unique |
| 19  | 2.75 x 10^11   | ~0.011            | Highly unique |
| 21  | 4.40 x 10^12   | ~7 x 10^-4        | Essentially unique |
| 25  | 1.13 x 10^15   | ~2.7 x 10^-6      | Unique (barring repeats) |
| 31  | 4.61 x 10^18   | ~6.7 x 10^-10     | Unique even in repeats |
| 37  | 1.88 x 10^22   | ~1.6 x 10^-13     | Extreme uniqueness |
| 43  | 7.70 x 10^25   | ~4 x 10^-17       | Redundantly unique |

These calculations apply to non-repetitive sequence. In practice, the human genome is approximately 50% repetitive (SINEs, LINEs, segmental duplications, tandem repeats). For repetitive regions, uniqueness requires longer k values -- often k >= 31 to avoid ambiguity in Alu elements (~300 bp, present at ~1.1 million copies) and other common repeats.

The KITSUNE tool (Pornputtapong et al., 2020) empirically determines optimal k-mer lengths using cumulative relative entropy, average common features, and observed common features. Their analysis found optimal k varies by organism complexity: ~11 for viruses, ~17 for bacteria, and ~34 for fungi. For human genome applications, k=31 is a well-established choice that balances uniqueness against computational cost.

### Sensitivity vs Specificity Tradeoff

**Specificity** increases monotonically with k because longer k-mers are more unique. A k-mer that maps to a single locus in the genome provides unambiguous signal for that locus. At k=31, essentially all k-mers from non-repetitive regions of the human genome are unique.

**Sensitivity** has a more complex relationship with k, driven by multiple interacting factors:

1. **Number of k-mers per read**: A read of length L produces (L - k + 1) k-mers. Longer k means fewer k-mers per read.
2. **Coverage per k-mer**: Each k-mer overlaps more reads when k is small, yielding higher per-k-mer counts.
3. **Variant signal strength**: For SNVs, exactly k k-mers are disrupted, but each disrupted k-mer has fewer supporting reads at larger k.
4. **Edge effects**: K-mers near fragment ends are lost. This is proportionally worse at larger k.

The net effect: there exists an optimal k for any given combination of variant type, coverage, and fragment length. The theoretical optimum balances the number of informative k-mers (favoring larger k for more disrupted k-mers per variant) against per-k-mer coverage (favoring smaller k).

### Diminishing Returns Above k=31

The marginal gain in specificity from increasing k beyond 31 is minimal for non-repetitive human genome regions. At k=31, the probability of a random match is already ~10^-10. Going to k=41 improves this to ~10^-16, but the practical impact is negligible since k=31 k-mers are already essentially unique outside of repeats.

Meanwhile, the cost of increasing k is real:
- Fewer k-mers per read (and per cfDNA fragment)
- Lower per-k-mer counts (critical at low VAF)
- More k-mers spanning variant junctions (problematic for indels where the variant disrupts more k-mers)

The exception is repetitive regions, where even k=31 may produce ambiguous k-mers. Here, k=41 or longer provides meaningful improvement. This motivates a variable-k approach.

## Impact on SNV Detection

### K-mers Disrupted by a Single Base Change

For a single nucleotide variant (SNV), the set of k-mers that differ between reference and variant alleles is exactly those k-mers whose window includes the variant position. If the variant is at genomic position p, then k-mers starting at positions (p - k + 1) through p all include position p. That gives exactly **k** disrupted k-mers.

```
Reference:  ...NNNNN[A]NNNNN...
                     ^
                     SNV position

Disrupted k-mers (k=5 example):
  Position p-4: [NNNN[A]]  (variant at rightmost position)
  Position p-3: [NNN[A]N]
  Position p-2: [NN[A]NN]
  Position p-1: [N[A]NNN]
  Position p:   [[A]NNNN]  (variant at leftmost position)

  = 5 disrupted k-mers (= k)
```

**Larger k = more disrupted k-mers = stronger signal** because more k-mers carry the variant base. However, each of those k-mers also needs sufficient coverage to be counted. The expected count of any variant k-mer at a given VAF is:

```
expected_count = depth * VAF * (1 - (k-1) / read_length)
```

The `(1 - (k-1)/read_length)` factor accounts for the fact that not all reads covering the variant position produce all k disrupted k-mers -- reads that start or end near the variant position produce fewer k-mers.

For concrete numbers with paired-end 150bp reads at 5000x raw coverage:

| k   | Disrupted k-mers | Expected count per k-mer at 0.1% VAF | Total variant signal |
|-----|------------------|---------------------------------------|----------------------|
| 21  | 21               | 5000 * 0.001 * (1 - 20/150) = 4.33   | 21 * 4.33 = 90.9    |
| 25  | 25               | 5000 * 0.001 * (1 - 24/150) = 4.20   | 25 * 4.20 = 105.0   |
| 31  | 31               | 5000 * 0.001 * (1 - 30/150) = 4.00   | 31 * 4.00 = 124.0   |
| 37  | 37               | 5000 * 0.001 * (1 - 36/150) = 3.80   | 37 * 3.80 = 140.6   |
| 43  | 43               | 5000 * 0.001 * (1 - 42/150) = 3.60   | 43 * 3.60 = 154.8   |

At 0.1% VAF, even k=43 yields expected counts of 3.6 per k-mer, which is above typical detection thresholds (count >= 2). The total variant signal (sum of counts across disrupted k-mers) increases with k because the linear growth in disrupted k-mers outweighs the sublinear loss in per-k-mer count. This favors larger k for SNV detection -- but only when per-k-mer counts remain above the minimum threshold.

At 0.05% VAF the picture shifts:

| k   | Expected count per k-mer at 0.05% VAF | Above threshold (>=2)? |
|-----|----------------------------------------|------------------------|
| 21  | 2.17                                   | Yes (marginal)         |
| 31  | 2.00                                   | Barely                 |
| 43  | 1.80                                   | No                     |

At the very lowest VAFs, smaller k becomes necessary to keep per-k-mer counts above the detection threshold.

## Impact on INDEL Detection

### K-mers Disrupted by Insertions and Deletions

For an insertion of n bases, the number of disrupted k-mers is (k + n), because the insertion affects all k-mers whose window spans the insertion breakpoint, plus the insertion introduces n new k-mers not present in the reference. However, the critical constraint is that the variant k-mers must be fully spanned by sequencing reads.

For a k-mer to span an insertion junction, the read must extend at least (k - offset) bases into the insertion from one side and at least offset bases on the other side. The effective coverage of junction-spanning k-mers drops with larger k:

```
For k=31, insertion of 15bp:
  - A k-mer centered on the junction needs 16 ref bases + 15 insert bases = 31 total
  - Only reads that span at least 31 bp across the junction contribute
  - With 150bp reads, most reads covering the region will span it

For k=31, insertion of 25bp:
  - Junction k-mers need up to 31 bases from the insert side alone
  - Only 6 of the 56 disrupted k-mers are fully within the insertion
  - Many junction k-mers require reads to span >31bp across the breakpoint

For k=31, insertion of 50bp:
  - Some k-mers are entirely within the insertion (have no reference anchor)
  - Others span the junction but require very specific read positioning
  - Coverage of individual junction k-mers is severely reduced
```

The thesis (Section 6.5.3) noted that insertions >15 bp (approaching k/2 for k=31) had significantly reduced k-mer coverage at the variant junction. This is because as the insertion grows relative to k, the fraction of disrupted k-mers that require reads spanning the full junction decreases.

**Shorter k dramatically helps INDEL detection**:

| k   | 10bp insertion | 20bp insertion | 30bp insertion |
|-----|----------------|----------------|----------------|
| 21  | 31 disrupted k-mers, all well-covered | 41 disrupted, most covered | 51 disrupted, junction k-mers still covered |
| 31  | 41 disrupted, all covered | 51 disrupted, junction coverage reduced | 61 disrupted, many junction k-mers poorly covered |
| 43  | 53 disrupted, some junction loss | 63 disrupted, significant junction loss | 73 disrupted, severe junction coverage loss |

For deletions, the situation is somewhat different. A deletion of n bases disrupts (k + n) k-mers in the reference, but creates only k junction k-mers in the variant allele. These junction k-mers are well-covered as long as reads span the deletion breakpoint, which is easier than spanning an insertion because the variant allele is shorter.

However, deletions in homopolymer contexts create ambiguous k-mers because multiple deletion positions produce the same junction k-mer sequence, confounding the walking algorithm.

### Thesis Recommendation: Test k=21, 25, 31, 37, 43

The thesis (Section 6.5.3) recommended testing five k values spanning the range from shorter (higher INDEL sensitivity) to longer (higher specificity). The odd values are conventional in k-mer analysis because they avoid palindromic k-mers (a k-mer cannot equal its own reverse complement when k is odd, simplifying canonical form computation).

## cfDNA Fragment Size Interaction

### Fragment Size Distribution

Cell-free DNA (cfDNA) from blood exhibits a characteristic fragment size distribution centered around nucleosomal wrapping:
- **Primary peak**: ~167 bp (mono-nucleosome + linker)
- **Secondary peak**: ~330 bp (di-nucleosome)
- **ctDNA enrichment**: Tumor-derived cfDNA fragments are enriched at shorter sizes (90-150 bp), with studies showing a 2-fold enrichment when selecting fragments in this range (Mouliere et al., 2018)

The fragment size directly constrains the number of k-mers extractable per fragment:

```
k-mers per fragment = fragment_length - k + 1
```

| Fragment length | k=21 | k=25 | k=31 | k=37 | k=43 |
|-----------------|------|------|------|------|------|
| 90 bp (short ctDNA) | 70 | 66 | 60 | 54 | 48 |
| 120 bp          | 100 | 96 | 90 | 84 | 78 |
| 150 bp          | 130 | 126 | 120 | 114 | 108 |
| 167 bp (mono-nucleosome) | 147 | 143 | 137 | 131 | 125 |
| 330 bp (di-nucleosome) | 310 | 306 | 300 | 294 | 288 |

For the primary cfDNA peak at 167 bp:
- k=21: 147 k-mers (88% of fragment used)
- k=31: 137 k-mers (82% of fragment used)
- k=43: 125 k-mers (75% of fragment used)

The proportional loss from k=21 to k=43 is (147-125)/147 = 15%. This is modest for the dominant 167 bp peak, but becomes more significant for the shorter ctDNA-enriched fragments:

For 90 bp fragments (enriched for ctDNA):
- k=21: 70 k-mers
- k=43: 48 k-mers (31% fewer)

Since ctDNA fragments tend to be shorter than non-tumor cfDNA, and since detecting ctDNA at low VAF is the primary goal, **larger k disproportionately impacts the signal from the molecules we most want to detect**.

### Edge Effects and Terminal K-mers

K-mers near fragment ends are particularly important because they may be the only k-mers covering a variant near the fragment boundary. With larger k, more of the fragment is "wasted" on terminal positions that cannot contribute complete k-mers. The first and last (k-1) bases of each fragment contribute to fewer than k overlapping k-mers.

For a 167 bp fragment, the fraction of bases contributing to the full k overlapping k-mers:
- k=21: (167 - 2*20) / 167 = 76%
- k=31: (167 - 2*30) / 167 = 64%
- k=43: (167 - 2*42) / 167 = 50%

At k=43, half the fragment's bases are in the "edge zone" where they contribute to fewer than k k-mers. This further reduces effective coverage for variant k-mers near fragment ends.

## Literature on Multi-k Approaches in Genome Assembly

### SPAdes Multi-k Strategy

SPAdes (Bankevich et al., 2012) is a de Bruijn graph assembler that uses multiple k-mer lengths iteratively. The approach works in stages:

1. Build a de Bruijn graph at a small k (e.g., k=21)
2. Use the graph to identify reliable paths and correct errors
3. Rebuild at a larger k (e.g., k=33), using information from the previous stage
4. Repeat at even larger k values (e.g., k=55, k=77)

The rationale is that smaller k values are better for:
- Low-coverage regions (more reads contribute k-mers)
- Error correction (errors create isolated tips that are easily identified)
- Gap filling (k-mers bridge regions where larger k would have insufficient coverage)

While larger k values are better for:
- Resolving repeats (longer k-mers are more unique)
- Reducing graph complexity (fewer spurious edges)
- Producing longer contigs

SPAdes auto-selects k values based on read length. For 150 bp reads, the default series is typically k=21,33,55,77.

### MEGAHIT Iterative Multi-k

MEGAHIT (Li et al., 2015) uses a strategy derived from IDBA assemblers where it iteratively builds succinct de Bruijn graphs (SdBGs) from small k to large k. The default k series is 21, 41, 61, 81, 99. At each k, contigs from the previous k iteration are used to guide graph construction, with small k favoring low-coverage regions and large k resolving repeats.

The key difference from SPAdes is that MEGAHIT uses succinct representations (compressed bit vectors) for the de Bruijn graph, enabling much lower memory usage at the cost of some speed.

### SKESA Strategic K-mer Extension

SKESA (Souvorov et al., 2018) takes a different approach: rather than building full de Bruijn graphs at multiple k values, it strategically extends contigs using different k values as needed. It starts with a hash-based de Bruijn graph at a moderate k, then extends contigs using longer k-mers to resolve ambiguities. This "extend as needed" philosophy is closer to what kmerdet could adopt -- using shorter k in regions where the standard k=31 fails to extend.

### Relevance to Variant Detection

The assembly literature provides two key insights for variant detection:

1. **Small k catches signal in low-coverage regions**: For low-VAF variants, the variant allele effectively has "low coverage." Using smaller k increases the number of reads contributing k-mers at the variant site, analogous to how assemblers use small k to bridge low-coverage gaps.

2. **Large k resolves ambiguity in repetitive contexts**: Variant targets in repetitive regions need longer k to ensure anchor k-mers are unique. This is directly analogous to repeat resolution in assembly.

3. **Iterative refinement**: Rather than running independently at each k, information from one k can inform analysis at another k. For variant detection, this could mean using k=21 results to identify candidate variant regions, then confirming with k=31 for specificity.

## Jellyfish Multi-k Considerations

### Running Multiple JF Databases

Jellyfish databases are tied to a single k value at creation time. To analyze at multiple k values, separate databases must be created:

```bash
jellyfish count -m 21 -s 100M -t 4 -C reads.fastq -o counts_k21.jf
jellyfish count -m 31 -s 100M -t 4 -C reads.fastq -o counts_k31.jf
jellyfish count -m 41 -s 100M -t 4 -C reads.fastq -o counts_k41.jf
```

### Resource Cost of Multi-k

**Time**: Each jellyfish count pass reads the entire FASTQ file. For 3 k values, counting time triples. With a typical liquid biopsy panel generating ~20M read pairs, each counting pass takes ~1.5 minutes, so 3 passes add ~3 minutes.

**Memory**: Each database is held independently. Memory per database is proportional to the number of distinct k-mers, which increases slightly with smaller k (more shared k-mers at smaller k). Rough estimates for targeted panel data after HUMID dedup:

| k   | Approximate distinct k-mers | Memory (at 16 bytes/entry) |
|-----|----------------------------|---------------------------|
| 21  | ~40M                       | ~640 MB                   |
| 31  | ~45M                       | ~720 MB                   |
| 41  | ~48M                       | ~768 MB                   |

For targeted panels (not WGS), multi-k is feasible on standard workstations with 8-16 GB RAM.

**Disk**: .jf files for targeted panel data are typically 200-800 MB each. Three databases require 0.6-2.4 GB of temporary disk space.

**Query time impact**: During k-mer walking, each extension step queries one database. With multi-k, if walking is performed independently at each k, query time scales linearly with the number of k values. However, each walk is independent and can be parallelized.

### Practical Considerations

For the kmerdet pipeline, multi-k analysis adds modest overhead to the counting phase but potentially significant overhead to the walking phase (multiple independent walks per target per k value). The key question is whether the improved detection sensitivity justifies the additional computation, which is addressed in the companion document on multi-k strategies.

## References

- Bankevich, A. et al. (2012). SPAdes: A new genome assembly algorithm. *Journal of Computational Biology*, 19(5), 455-477.
- Li, D. et al. (2015). MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics*, 31(10), 1674-1676.
- Souvorov, A. et al. (2018). SKESA: strategic k-mer extension for scrupulous assemblies. *Genome Biology*, 19, 153.
- Pornputtapong, N. et al. (2020). KITSUNE: A Tool for Identifying Empirically Optimal K-mer Length. *Frontiers in Bioengineering and Biotechnology*, 8, 556413.
- Mouliere, F. et al. (2018). Enhanced detection of circulating tumor DNA by fragment size analysis. *Science Translational Medicine*, 10(466).
- Marcais, G. & Kingsford, C. (2011). A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. *Bioinformatics*, 27(6), 764-770.
- Kokot, M. et al. (2017). KMC 3: counting and manipulating k-mer statistics. *Bioinformatics*, 33(17), 2759-2761.
