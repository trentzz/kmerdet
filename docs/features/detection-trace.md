# Detection Trace — Verbose Per-Target Analysis

## What It Does

Generates a detailed trace of the detection process for each target, showing every step: walking exploration, graph construction, path finding, classification, and quantification. Allows researchers to understand exactly why a variant was or wasn't found.

## Why It Matters

When a variant is missed or incorrectly classified, researchers need to diagnose where in the pipeline the signal was lost:
- Did the walker explore the variant k-mers? (walking issue)
- Were the variant k-mers connected in the graph? (graph construction issue)
- Was a path found through variant nodes? (pathfinding issue)
- Was the path correctly classified? (classification issue)
- Was the rVAF too low to report? (quantification issue)

This is inherently verbose output — it's a debugging/analysis tool, not standard output.

## CLI Interface

```bash
# Generate trace for all targets
kmerdet detect --db sample.jf --targets targets/ --trace trace_output/

# Generate trace for specific targets only
kmerdet detect --db sample.jf --targets targets/ --trace trace_output/ --trace-targets "TP53_R248W,KRAS_G12D"

# Trace as part of run pipeline
kmerdet run --db sample.jf --targets-dir targets/ --expected expected.tsv --trace trace_output/
```

## Output Structure

```
trace_output/
  summary.json              # Overview of all targets
  targets/
    TP53_R248W.json          # Full trace for this target
    KRAS_G12D.json
    ...
```

## Trace Contents Per Target

### Walking Phase
```json
{
  "walking": {
    "reference_kmers": 35,
    "reference_kmer_counts": [1205, 1198, ...],
    "total_nodes_discovered": 42,
    "alt_nodes_discovered": 7,
    "branching_points": 2,
    "max_stack_depth_reached": 12,
    "limits_hit": null,  // or "max_node" / "max_stack" / "max_break"
    "extensions": [
      {
        "from": "ACGTACGTACGTACGTACGTACGTACGTACG",
        "children": [
          {"sequence": "CGTACGTACGTACGTACGTACGTACGTACGT", "count": 1200, "is_reference": true},
          {"sequence": "CGTACGTACGTACGTACGTACGTACGTACGA", "count": 15, "is_reference": false}
        ],
        "branching": true
      }
    ]
  }
}
```

### Graph Phase
```json
{
  "graph": {
    "total_nodes": 44,
    "reference_nodes": 35,
    "alt_nodes": 7,
    "virtual_nodes": 2,
    "total_edges": 52,
    "reference_edges": 34,
    "alt_edges": 18
  }
}
```

### Pathfinding Phase
```json
{
  "pathfinding": {
    "paths_found": 2,
    "reference_path": {
      "length": 35,
      "sequence": "ACGT...",
      "total_weight": 0.34
    },
    "alternative_paths": [
      {
        "index": 1,
        "length": 36,
        "sequence": "ACGT...T...",
        "total_weight": 3.21,
        "divergence_point": 15,
        "rejoin_point": 22,
        "shared_prefix_length": 15,
        "shared_suffix_length": 13
      }
    ]
  }
}
```

### Classification Phase
```json
{
  "classification": [
    {
      "path_index": 1,
      "variant_type": "Insertion",
      "variant_name": "15:A/ACGT:19",
      "ref_allele": "A",
      "alt_allele": "ACGT",
      "start": 15,
      "end": 19,
      "diff_details": {
        "left_trim": 15,
        "right_trim": 13,
        "ref_length_after_trim": 1,
        "alt_length_after_trim": 4
      }
    }
  ]
}
```

### Quantification Phase
```json
{
  "quantification": {
    "num_paths": 2,
    "num_unique_kmers": 42,
    "contribution_matrix_shape": [42, 2],
    "coefficients": [1180.5, 14.2],
    "rvafs": [0.988, 0.012],
    "min_coverages": [1190, 12],
    "nnls_iterations": 1,
    "residual_norm": 45.2
  }
}
```

## Data Structures

```rust
pub struct DetectionTrace {
    pub target: String,
    pub walking: WalkingTrace,
    pub graph: GraphTrace,
    pub pathfinding: PathfindingTrace,
    pub classifications: Vec<ClassificationTrace>,
    pub quantification: Option<QuantificationTrace>,
    pub outcome: String,  // "variant_detected" / "reference_only" / "error"
}

pub struct WalkingTrace {
    pub reference_kmers: usize,
    pub reference_kmer_counts: Vec<u64>,
    pub total_nodes: usize,
    pub alt_nodes: usize,
    pub branching_points: usize,
    pub max_stack_depth: usize,
    pub limits_hit: Option<String>,
    pub extensions: Vec<ExtensionTrace>,
}

pub struct ExtensionTrace {
    pub from: String,
    pub children: Vec<ChildTrace>,
    pub branching: bool,
}
```

## Acceptance Criteria

- [ ] `--trace <dir>` flag on detect and run subcommands
- [ ] Per-target JSON trace files in output directory
- [ ] summary.json with overview of all targets
- [ ] Walking trace: nodes discovered, branching, limits hit
- [ ] Graph trace: node/edge counts by type
- [ ] Pathfinding trace: all paths with divergence/rejoin points
- [ ] Classification trace: diff details for each variant
- [ ] Quantification trace: NNLS details, coefficients, rVAFs
- [ ] `--trace-targets` filter to limit trace to specific targets
- [ ] Extension-level detail (which k-mers were explored and their counts)
- [ ] Trace output does not impact detection performance when not enabled
