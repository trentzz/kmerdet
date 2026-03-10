# km Tool Analysis (iric-soft/km)

## Overview

km is a Python-based bioinformatics tool for **targeted variant detection using k-mer decomposition** of RNA-seq (or DNA-seq) data. Rather than performing genome-wide variant calling, km focuses on specific regions of interest defined by user-provided target sequences.

- **Language**: Python (95.8%), Shell (4.2%)
- **License**: MIT
- **PyPI Package**: `km-walk`
- **Repository**: https://github.com/iric-soft/km
- **Python Support**: 3.8-3.12
- **Jellyfish Requirement**: 2.2+ with Python bindings

## What It Detects

- Single nucleotide variants (SNVs)
- Insertions
- Deletions
- Internal tandem duplications (ITDs)
- Complex indels and multi-nucleotide variants (MNVs)

## Architecture & Module Structure

```
km/
├── __init__.py
├── __main__.py          # CLI entry point
├── km.py                # Core module
├── argparser/
│   ├── common.py
│   ├── find_mutation.py
│   ├── find_report.py
│   ├── linear_kmin.py
│   └── min_cov.py
├── tools/
│   ├── find_mutation.py  # Primary analysis tool
│   ├── find_report.py    # Post-processing/reporting
│   ├── linear_kmin.py    # Minimum k-length calculator
│   └── min_cov.py        # Coverage statistics
└── utils/
    ├── Graph.py           # Directed weighted graph + Dijkstra-like pathfinding
    ├── Jellyfish.py       # Python interface to jellyfish DB
    ├── MutationFinder.py  # Core algorithm: k-mer walking + variant identification
    ├── PathQuant.py       # Path quantification via linear regression
    ├── Sequence.py        # Reference/alternative sequence objects
    └── common.py          # Shared utilities
```

### Pre-built Target Catalogs

```
data/catalog/
├── GRCh37/   (hg19)
│   ├── DNMT3A_R882_exon_23.fa
│   ├── FLT3-ITD_exons_13-15.fa
│   ├── FLT3-TKD_exon_20.fa
│   ├── IDH1_R132.fa
│   ├── KMT2A-PTD_8-2.fa
│   ├── NPM1_4ins_exons_10-11utr.fa
│   ├── NSD1_exon6-NUP98_exon13.fa
│   └── NUP98_exon11-NSD1_exon7.fa
└── GRCh38/   (hg38 - same targets)
```

These are all **leukemia-associated mutations** (AML/MDS), which is km's primary application domain.

## Core Algorithm: How find_mutation Works

### Step 1: Initialization (MutationFinder.__init__)

1. Load the target FASTA sequence and decompose it into k-mers (`RefSeq`)
2. Open the jellyfish database (`Jellyfish`)
3. Query the jellyfish DB for **every k-mer in the reference sequence**
4. For each reference k-mer, perform **depth-first k-mer extension** (`__extend`)

### Step 2: K-mer Extension (DFS Walking)

The `__extend` method performs recursive depth-first search from each reference k-mer:

```python
def __extend(self, stack, breaks=0):
    cur_seq = stack[-1]
    childs = self.jf.get_child(cur_seq, forward=True)  # Get all 4 possible extensions

    for child in childs:
        if child not in node_data and child not in stack:
            self.__extend(stack + [child], breaks)  # Recurse
```

**Child extension logic** (`Jellyfish.get_child`):
- For a k-mer of length k, try appending each of {A, C, G, T} and taking the last k bases
- Query jellyfish for the count of each candidate
- Keep only candidates whose count >= `max(sum * cutoff, n_cutoff)`
  - Default `cutoff=0.30` (30% of total sibling counts)
  - Default `n_cutoff=500` (absolute minimum count)
- This thresholding prevents chasing sequencing errors

**Limits**:
- `max_stack=500`: maximum depth of DFS extension
- `max_break=10`: maximum number of branching points
- `max_node=10000`: maximum total nodes

### Step 3: Graph Construction (graph_analysis)

1. Build a directed weighted graph where each node is a k-mer
2. Edges connect k-mers that overlap by (k-1) bases
3. **Reference edges get weight 0.01**, non-reference edges get weight 1.0
4. Add virtual "BigBang" (source) and "BigCrunch" (sink) nodes
5. Use Dijkstra-like shortest path algorithm from both directions
6. **Remove reference edges** from the graph
7. Find **all shortest paths** through remaining non-reference edges

This produces all alternative paths (variant sequences) that differ from the reference.

### Step 4: Variant Identification (diff_path_without_overlap)

Compare reference path to each alternative path:
- Walk from the left to find the first differing k-mer position (`i`)
- Walk from the right to find the last differing k-mer position (`j_ref`, `j_seq`)
- Classify the variant:
  - `start == end_ref == end_var` → **Reference** (no change)
  - `end_ref == end_var` → **Substitution** (SNV/MNP)
  - `start == end_ref_overlap` → **ITD** (internal tandem duplication)
  - `end_ref < end_var` with no deleted bases → **Insertion**
  - `end_ref > end_var` with no inserted bases → **Deletion**
  - Otherwise → **Indel** (complex)

### Step 5: Quantification (PathQuant)

Uses **linear regression** to estimate the expression level of each path:

1. Build a contribution matrix: `contrib[kmer_i, path_j]` = number of times k-mer i appears in path j
2. Solve `contrib @ coef = counts` via least squares (`np.linalg.lstsq`)
3. Refine with gradient descent to eliminate negative coefficients
4. Compute **rVAF** (relative variant allele frequency): `coef / sum(coef)`

This elegantly handles overlapping k-mers between paths by decomposing the observed counts into per-path contributions.

### Step 6: Cluster Quantification

When multiple mutations overlap in k-mer space, they're grouped into clusters and quantified together. This handles cases where compound heterozygous mutations would cause zero minimum coverage for all individual paths.

## Input/Output

### Input
- **Target FASTA**: User-designed sequences (few hundred bp) around regions of interest
- **Jellyfish DB** (`.jf`): Pre-computed k-mer count database from sequencing reads

### Output (TSV columns)
| Column | Description |
|--------|-------------|
| Database | Jellyfish DB filename |
| Query | Target sequence name |
| Type | Variant type (Reference/Substitution/Insertion/Deletion/ITD/Indel) |
| Variant_name | Position and change (e.g., `41:acgt/ACGTACGT:45`) |
| rVAF | Relative variant allele frequency (0.0-1.0) |
| Expression | Estimated expression coefficient |
| Min_coverage | Minimum k-mer count along the path |
| Start_offset | Offset for cluster mode |
| Sequence | Full variant sequence |
| Reference_expression | Reference path expression |
| Reference_sequence | Full reference sequence |
| Info | Additional notes (vs_ref or cluster info) |

## Key Parameters

- **k-mer length**: Minimum 21, default 31 (higher = more specific, lower sensitivity near edges)
- **ratio** (cutoff): 0.30 — minimum fraction of sibling k-mer counts for extension
- **count** (n_cutoff): 500 — absolute minimum k-mer count for extension
- **steps** (max_stack): 500 — maximum DFS depth
- **branchs** (max_break): 10 — maximum branching points in DFS

## Key Design Decisions

1. **Targeted, not genome-wide**: Requires pre-designed target sequences → fast, focused
2. **Graph-based pathfinding**: Weighted graph with Dijkstra ensures shortest/most likely paths
3. **Linear regression for VAF**: Elegant decomposition of k-mer counts into per-allele contributions
4. **Threshold-based extension**: Prevents chasing sequencing errors while allowing low-VAF detection
5. **Cluster mode**: Handles compound mutations that would confuse simple path quantification

## References

- Audoux et al. (2017) — Original km paper (leukemia variant detection)
- Audoux et al. (2019) — Extended km methodology
- Focus on AML: NPM1, FLT3-ITD, DNMT3A, IDH1/2 mutations
