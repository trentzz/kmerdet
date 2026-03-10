# K-mer Variant Detection Algorithm

## Fundamental Concept

Instead of aligning reads to a reference genome, k-mer variant detection works by:

1. Pre-counting all k-mers in sequencing data (via jellyfish)
2. Defining a target region around a suspected mutation
3. Decomposing the target into overlapping k-mers
4. Walking through k-mer space to find all paths supported by the data
5. Comparing alternative paths to the reference to identify variants
6. Quantifying variant allele frequency from k-mer count decomposition

## The Algorithm in Detail

### Phase 1: K-mer Database Preparation

```
Sequencing reads → jellyfish count → .jf database
```

The jellyfish database maps every k-mer (substring of length k) in the sequencing data to its occurrence count. For a 31-mer database from a typical WGS sample:
- ~3 billion k-mers extracted from reads
- ~2-3 billion distinct k-mers stored
- Count reflects sequencing depth at that k-mer

### Phase 2: Target Sequence Design

A target sequence is a short (typically 150-500 bp) reference sequence surrounding a known or suspected mutation site. Requirements:

- Must be long enough that the first k-mer and last k-mer are **unique** in the genome
- The first and last k-mers serve as **anchor points** (source and sink)
- All possible variant sequences must be constructable between these anchors

Example: For NPM1 4bp insertion detection, the target covers ~200bp around exon 10-11.

### Phase 3: K-mer Walking (Graph Construction)

Starting from each reference k-mer, perform depth-first extension:

```
Reference k-mers:  [ACGT...31bp...] → [CGT...31bp...A] → [GT...31bp...AC] → ...
                                    ↘ [CGT...31bp...T]    (branch = variant!)
```

**Extension rule**: For k-mer `S[1..k]`, try all 4 extensions `S[2..k]+{A,C,G,T}`. Keep any extension whose count in the jellyfish DB exceeds the threshold.

**Thresholding**:
```
threshold = max(sum_of_sibling_counts * 0.30, absolute_minimum)
```
- This relative threshold adapts to local coverage
- The absolute minimum prevents following noise at high-coverage loci

### Phase 4: Graph Pathfinding

The k-mer graph is a directed weighted graph:
- **Nodes**: All k-mers discovered during walking
- **Edges**: k-mer A → k-mer B if A[2..k] == B[1..k-1] (k-1 overlap)
- **Weights**: Reference edges = 0.01, non-reference edges = 1.0

Algorithm:
1. Run Dijkstra from source ("BigBang") → find predecessor for each node
2. Run Dijkstra from sink ("BigCrunch") → find successor for each node
3. Remove all reference edges
4. For each remaining non-reference edge, construct the shortest path through it
5. Deduplicate paths → these are the **alternative (variant) paths**

### Phase 5: Variant Classification

Compare reference path `R` to alternative path `A`:

```
R: [k1][k2][k3][k4][k5][k6][k7][k8][k9][k10]
A: [k1][k2][k3][kX][kY][k8][k9][k10]
                  ^--- divergence ---^
```

Walk from left: find first position where R[i] != A[i] → `start`
Walk from right: find last position where they differ → `end_ref`, `end_var`

Classification logic:
| Condition | Variant Type |
|-----------|-------------|
| `start == end_ref == end_var` | Reference (identical) |
| `end_ref == end_var` (same length) | Substitution (SNV/MNP) |
| `start == end_ref_overlap` | Internal Tandem Duplication (ITD) |
| `end_ref < end_var`, no deleted bases | Insertion |
| `end_ref > end_var`, no inserted bases | Deletion |
| Otherwise | Complex Indel |

### Phase 6: VAF Quantification

The key insight: each k-mer's count is the **sum** of contributions from all paths passing through it.

**Mathematical model**:
```
For each k-mer i:  count[i] = Σ(contrib[i,j] * coef[j])  for all paths j
```

Where:
- `count[i]` = observed jellyfish count for k-mer i
- `contrib[i,j]` = 1 if k-mer i appears in path j, 0 otherwise (or >1 for ITDs)
- `coef[j]` = expression level (coverage) of path j

**Solution**: Solve via least squares regression:
```
coef = argmin ||contrib @ coef - counts||²
```

Then refine with gradient descent to eliminate negative coefficients.

**rVAF** (relative variant allele frequency):
```
rVAF[j] = coef[j] / Σ(coef)
```

This gives the fraction of reads supporting each path, directly analogous to traditional VAF.

## Handling Specific Variant Types

### SNVs (Single Nucleotide Variants)
- Exactly `k` k-mers will differ between reference and variant paths
- The mutated base is at the center of these differing k-mers
- Both paths have the same length

### Insertions
- The variant path is longer than the reference path
- The number of additional k-mers equals the insertion length
- New k-mers span the insertion breakpoints

### Deletions
- The variant path is shorter than the reference path
- Some reference k-mers are missing from the variant path
- New junction k-mers span the deletion breakpoints

### Internal Tandem Duplications (ITDs)
- The variant path contains a duplicated segment
- Some k-mers appear twice in the variant path (`contrib > 1`)
- This is handled by the contribution matrix allowing values > 1
- Detected when `start == end_ref_overlap` (the whole reference is retraced)

### Complex/Compound Mutations
- When multiple mutations overlap in k-mer space, they form a "cluster"
- Cluster mode clips the sequence around each mutation group
- Quantifies all paths in the cluster simultaneously
- Prevents mutual interference between overlapping mutations

## Advantages Over Alignment-Based Approaches

1. **No reference bias**: K-mers directly represent the sequencing data, not alignment interpretations
2. **Better for complex variants**: ITDs, large indels, and structural variants that are hard to align
3. **Speed**: Pre-computed k-mer database enables fast targeted queries
4. **Sensitivity at low VAF**: Direct k-mer counting can detect variants at lower frequencies than alignment-based callers
5. **No mapping quality issues**: Avoids problems with multi-mapped reads and mapping quality filtering
6. **Handles repetitive regions**: K-mers spanning variant junctions are kept regardless of mapping status

## Limitations

1. **Targeted only**: Requires pre-designed target sequences (not genome-wide discovery)
2. **K-mer length constraint**: Insertions/deletions larger than k/2 have reduced sensitivity
3. **Requires sufficient coverage**: Low-count k-mers are filtered out as noise
4. **Homopolymer regions**: K-mers in repetitive regions may not be unique
5. **Adjacent mutations**: Multiple nearby mutations can complicate path finding

## Key Parameters and Their Effects

| Parameter | Default | Effect |
|-----------|---------|--------|
| k (k-mer length) | 31 | Higher = more specific, lower sensitivity for large indels |
| ratio (cutoff) | 0.30 | Lower = more sensitive but more noise, higher = more specific |
| count (n_cutoff) | 500 | Lower = detect lower VAF, higher = reduce false positives |
| max_stack | 500 | Maximum depth of k-mer extension |
| max_break | 10 | Maximum branching points (limits search space) |

## References

- Audoux et al. — km: RNA-seq investigation using k-mer decomposition
- KmerVC (Stanford) — K-mer validation of cancer mutations
- 2-kupl — Mapping-free variant detection from DNA-seq
- DE-kupl — Exhaustive capture of biological variation via k-mer decomposition
