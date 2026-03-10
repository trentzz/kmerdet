# Graph Construction and Pathfinding Improvements

## 1. Current Weighting Scheme Analysis

### The Binary Weight Model

The current graph construction in `src/graph/builder.rs` uses two fixed edge weights:

```rust
pub const REF_EDGE_WEIGHT: f64 = 0.01;
pub const ALT_EDGE_WEIGHT: f64 = 1.0;
```

All reference edges get weight 0.01, all non-reference edges get weight 1.0. This creates a 100x preference for reference paths in Dijkstra's algorithm. The consequence: Dijkstra always finds the reference path as the shortest path, which is by design -- the algorithm then removes reference edges and searches for shortest alternative paths.

### Why These Specific Values?

The 0.01/1.0 values appear empirical rather than theoretically motivated. From the km source code, they were likely chosen to ensure:

1. The reference path is always the globally shortest path (any sequence of reference edges is shorter than any single non-reference edge detour)
2. Among alternative paths, Dijkstra finds the shortest by hop count (all non-reference edges are equal weight)

This works for the core use case but has limitations:

**No biological signal in weights.** A non-reference edge supported by 5,000 reads and one supported by 2 reads have the same weight. Dijkstra cannot distinguish well-supported variant paths from noise paths.

**No path ranking by evidence.** When multiple alternative paths exist, they are all reported regardless of support. The user must manually examine rVAF values to determine which are real.

**Assumption of single best path.** The binary scheme assumes we want the reference path plus deviations from it. It cannot naturally handle regions where two non-reference alleles are both common (e.g., triallelic sites).

## 2. Alternative Weighting Approaches

### Coverage-Weighted Edges

Replace fixed weights with coverage-derived weights:

```
weight(A -> B) = 1 / count(B)
```

Higher count means lower weight, so Dijkstra prefers high-coverage paths. The reference path naturally gets the lowest weight because reference k-mers typically have the highest counts.

**Advantages:**
- Paths are ranked by coverage support without requiring a separate quantification step
- Variant paths at reasonable VAF get moderate weight, distinguishing them from noise
- Naturally handles varying coverage across the target

**Disadvantages:**
- Log-scale counts mean the weight difference between count=5000 and count=50 is only 100x, while the biological significance difference is much larger
- Count=1 k-mers get weight=1.0, which may be comparable to genuine low-VAF variants

**Practical concern:** The reference/non-reference edge distinction is still needed for the algorithm's step 3 (remove reference edges, find alternative paths). Coverage weighting should augment, not replace, the reference classification.

### Log-Coverage Weighting

```
weight(A -> B) = -log2(count(B) / total_coverage)
```

This has an information-theoretic interpretation: weight is the surprisal of observing that k-mer. A k-mer at expected coverage has low surprisal (low weight), while a rare k-mer has high surprisal (high weight).

**Example at 5000x coverage:**

| K-mer type | Count | Weight |
|-----------|-------|--------|
| Reference | 5000 | -log2(5000/5000) = 0.0 |
| Variant at 10% VAF | 500 | -log2(500/5000) = 3.32 |
| Variant at 1% VAF | 50 | -log2(50/5000) = 6.64 |
| Variant at 0.1% VAF | 5 | -log2(5/5000) = 9.97 |
| Error | 2 | -log2(2/5000) = 11.29 |

The log scale compresses the range while maintaining rank order. The gap between a 0.1% VAF variant and an error k-mer is only ~1.3, which may be insufficient for reliable discrimination. This approach works best when combined with a minimum count threshold to pre-filter noise.

### Ratio-Based Weighting

```
weight(A -> B) = 1 - min(1, count(B) / expected_count)
```

Where `expected_count` is the local median coverage from reference k-mers. K-mers at expected coverage get weight ~0, k-mers far below expected coverage get weight ~1.

This normalizes for coverage variation across the target. A k-mer with count=500 in a 5000x region and a k-mer with count=50 in a 500x region both get the same weight (0.9), which is appropriate because both represent ~10% VAF.

### Recommended Hybrid Approach

Retain the binary reference/non-reference classification for the core algorithm (steps 1-3: Dijkstra, remove reference edges, find alternative paths). Then, in step 4, rank the discovered alternative paths by a coverage-based score:

```rust
fn path_score(path: &KmerPath, graph: &KmerGraph) -> f64 {
    let counts: Vec<u64> = path.kmers.iter()
        .map(|k| graph.nodes[graph.node_index[k]].count)
        .collect();
    // Use minimum count as the bottleneck metric
    *counts.iter().min().unwrap() as f64
}
```

This decouples pathfinding (structural) from ranking (quantitative), keeping the algorithm simple while providing meaningful path ordering.

## 3. Handling Complex and Compound Variants

### Multiple Variants Within K-mer Range

When two variants are within `k` bases of each other, their k-mer perturbations overlap. The k-mers affected by variant A include some that are also affected by variant B. This creates a complex graph topology with multiple branching and merging points.

**Example:** SNV at position 100 and SNV at position 120, with k=31. K-mers spanning positions 70-130 are affected by one or both variants. The graph has four possible paths through this region:

1. Reference + Reference (wild type)
2. Variant A + Reference
3. Reference + Variant B
4. Variant A + Variant B (compound)

The number of paths grows exponentially with the number of nearby variants: `n` variants produce up to `2^n` paths.

### Current Cluster Approach

km handles this via cluster mode: when multiple mutations overlap in k-mer space, they are grouped and quantified jointly. The contribution matrix includes all `2^n` paths (or those actually observed in the data), and NNLS decomposes the k-mer counts into per-path coefficients.

This works well for 2-3 nearby variants but becomes computationally expensive for larger clusters. The contribution matrix has dimensions `(number_of_kmers x number_of_paths)`, and with `2^n` paths the matrix grows exponentially.

### Proposed Improvements

**Early cluster detection during walking:** Instead of discovering all paths and clustering after, detect clusters during graph construction. When a branching point is found within `k` positions of a previous branching point, mark the region as a cluster and adjust walking parameters (e.g., increased max_break).

**Phased vs. unphased representation:** For liquid biopsy, variants on the same cfDNA fragment are in cis (same allele). Variants on different fragments may be in trans. K-mer walking cannot distinguish phasing directly, but the count patterns can:

- If paths AB and ab exist but not Ab or aB, the variants are phased (cis or trans on two alleles)
- If all four combinations exist, the variants may be independently segregating

Report phasing information in the output when detectable.

**Practical limit on cluster size:** Cap cluster enumeration at 2^8 = 256 paths. Beyond this, the NNLS decomposition is numerically unstable and the biological interpretation is ambiguous. Flag such regions for manual review.

## 4. Graph Pruning Strategies

Graph pruning removes noise before pathfinding, improving both performance and accuracy. The following strategies are adapted from de Bruijn graph assembly practices (Velvet, SPAdes, ABySS), modified for the variant detection context.

### Dead-End Removal

Paths that do not reach the sink node cannot represent complete variant sequences. After graph construction, perform reachability analysis:

1. BFS/DFS forward from source: mark all reachable nodes
2. BFS/DFS backward from sink: mark all reachable nodes
3. Remove any node not reachable from both source and sink

This is already implicit in the current algorithm (Dijkstra from both ends identifies reachable nodes), but explicit removal before pathfinding reduces the graph size and speeds up subsequent steps.

**Complexity:** O(V + E) for two BFS passes, negligible compared to pathfinding.

### Tip Clipping

Tips are short branches that extend from the graph but terminate without connecting to another path. In assembly, tips shorter than `2k` are removed. For variant detection, the threshold should be more conservative:

- Tips shorter than `k/2` extensions: remove (almost certainly errors)
- Tips between `k/2` and `k` extensions: flag but retain (could be truncated variant paths)
- Tips longer than `k` extensions: retain (possible novel sequence)

```rust
fn clip_tips(graph: &mut KmerGraph, k: usize) {
    let tip_threshold = k / 2;
    loop {
        let tips: Vec<usize> = graph.find_tips(tip_threshold);
        if tips.is_empty() { break; }
        for node_idx in tips {
            graph.remove_node(node_idx);
        }
    }
}
```

The loop is necessary because removing a tip may expose new tips (a branch that was connected only through the removed tip).

### Bubble Collapsing

Bubbles are pairs of short paths with the same start and end nodes. In the variant detection graph:

- **Error bubbles:** 1-2 base differences, one path has much lower coverage. Collapse to higher-coverage path.
- **Variant bubbles:** Diverge for `k` or more positions, both paths have appreciable coverage. Preserve.

The Velvet Tour Bus algorithm detects bubbles by running a modified Dijkstra from each node and identifying convergence points. For kmerdet's smaller graphs (typically <10,000 nodes), a simpler approach works:

```rust
fn find_bubbles(graph: &KmerGraph) -> Vec<Bubble> {
    let mut bubbles = Vec::new();
    for node in graph.branching_nodes() {
        for (path_a, path_b) in node.divergent_path_pairs() {
            if path_a.rejoins(&path_b) && path_a.len() < k {
                bubbles.push(Bubble { path_a, path_b });
            }
        }
    }
    bubbles
}
```

**Collapse criteria:** Remove the lower-coverage path only if:
1. Both paths are shorter than `k` k-mers (too short to be a real variant)
2. The lower-coverage path has less than 10% of the higher-coverage path's minimum count

### Low-Coverage Path Removal

After constructing the graph but before pathfinding, remove edges where the target node has count below a noise threshold. This is similar to the walking threshold but applied globally to the graph:

```
noise_threshold = max(2, median_coverage * error_rate)
```

Nodes with count below this threshold that are not reference k-mers are removed. This eliminates sporadic error k-mers that passed the walking threshold due to local count fluctuations.

### Pruning Order

Apply pruning in this order for best results:

1. **Dead-end removal** (removes disconnected nodes)
2. **Low-coverage path removal** (removes noise edges, may create new dead-ends)
3. **Dead-end removal again** (clean up after step 2)
4. **Tip clipping** (remove short dead-end branches)
5. **Bubble collapsing** (merge redundant error paths)

## 5. Alternative Pathfinding Algorithms

### Current: Dijkstra + All-Shortest-Paths

The current algorithm (from km):

1. Dijkstra forward from source ("BigBang") to all nodes
2. Dijkstra backward from sink ("BigCrunch") to all nodes
3. Remove all reference edges
4. For each non-reference edge (u, v), compute the shortest path through it: `dist_source_to_u + weight(u,v) + dist_v_to_sink`
5. Enumerate all paths with minimum total weight

This is efficient for finding all shortest paths but has limitations:
- Only finds paths tied for minimum weight (misses second-best paths)
- Cannot handle graphs with negative edges (not an issue with current weights)
- The "all shortest paths through non-reference edges" enumeration can be exponential in the number of non-reference edges

### A* Search

A* extends Dijkstra with a heuristic function `h(n)` that estimates the remaining distance from node `n` to the goal:

```
f(n) = g(n) + h(n)
```

Where `g(n)` is the known shortest distance from source to `n`. A* explores nodes in order of `f(n)`, prioritizing nodes that appear closer to the goal.

The AStarix algorithm (Ivanov et al., 2020) demonstrated that A* with domain-specific heuristics can be 1-2 orders of magnitude faster than standard shortest-path algorithms for sequence-to-graph alignment. Their heuristic uses the remaining query sequence to estimate the minimum edit distance to the goal.

**Heuristic for kmerdet:** For a node at position `p` in the reference, the minimum distance to the sink is approximately `(ref_length - p) * REF_EDGE_WEIGHT`, since the best case is traversing the remaining reference. This is admissible (never overestimates) and consistent (satisfies the triangle inequality), guaranteeing optimal results.

**When A* helps:** A* is most beneficial when the goal is to find a single shortest path quickly, not to enumerate all paths. For kmerdet's "find all alternative paths" objective, A* provides less benefit because we need exhaustive enumeration, not early termination.

**Recommendation:** Use A* for a potential "fast mode" that reports only the single best alternative path, useful for quick screening when full enumeration is not needed.

### K-Shortest Paths (Yen's Algorithm)

Yen's algorithm finds the top-K shortest paths between source and sink, ranked by total weight. It works by:

1. Find the shortest path P1 (Dijkstra)
2. For k = 2 to K:
   - For each node in P_{k-1}, temporarily remove that node's outgoing edge used by P_{k-1}
   - Find the shortest path in the modified graph
   - The best such path is P_k

**Complexity:** O(K * V * (V * log V + E)) -- K iterations of Dijkstra.

**Application to variant detection:** Instead of the binary "find all shortest paths through non-reference edges," find the top-K paths by weight. This provides a natural ranking of variant paths by support:

| Path rank | Interpretation |
|-----------|---------------|
| P1 | Reference path (shortest) |
| P2 | Most strongly supported variant |
| P3 | Second variant (or error path) |
| ... | Diminishing support |
| PK | Weakest reported path |

Choose K based on the expected number of alleles: K=3 for diploid samples, K=5-10 for mixed samples or liquid biopsy (where multiple subclones may be present).

**Advantage over current approach:** Naturally limits output to the most supported paths. The current approach can report hundreds of paths in repetitive regions, most of which are noise.

**Disadvantage:** Yen's algorithm finds simple paths (no repeated nodes). If a variant path must traverse a node used by the reference path (common for SNVs where surrounding k-mers are shared), Yen's may miss it. The current approach avoids this by removing reference edges entirely, creating a different graph topology.

**Recommendation:** Use Yen's algorithm as an optional alternative mode, particularly useful when graph pruning is insufficient to control path count.

### Eppstein's Algorithm

Eppstein's algorithm (1998) finds K shortest paths in O(E + V log V + K) time, dramatically faster than Yen's O(K * V * (V log V + E)). It allows paths with repeated vertices (not just simple paths), which is important for ITD detection where the variant path loops through reference nodes.

For kmerdet's typical graph sizes (V < 10,000, E < 40,000, K < 50), both algorithms are fast enough, but Eppstein's is preferable for its theoretical guarantees and support for non-simple paths.

### Flow-Based Approaches

**Max-flow/min-cut** can identify variant boundaries: the minimum cut between two reference nodes defines the set of edges that must be traversed by any path between them. If the min-cut passes through non-reference edges, those edges are part of a variant path.

**Application:** Finding the boundaries of an INDEL in the graph:

```
[reference k-mers]---min-cut---[INDEL region]---min-cut---[reference k-mers]
```

The min-cut provides a principled way to identify where variants start and end, complementing the sequence-based classification in `diff_path_without_overlap`.

**Practical concern:** Flow algorithms require edge capacities. The natural capacity is k-mer count, but count varies across the path. Use minimum count along each path as the capacity.

### Comparison Table

| Algorithm | Time Complexity | Output | Repeated Nodes | Best For |
|-----------|----------------|--------|----------------|----------|
| Dijkstra + all-shortest | O(V log V + E) + enumeration | All paths tied for shortest | N/A (ref edges removed) | Current approach, complete enumeration |
| A* | O(V log V + E) with heuristic | Single shortest path | No | Fast screening mode |
| Yen's K-shortest | O(K * V * (V log V + E)) | Top-K simple paths | No | Ranked variant reporting |
| Eppstein's K-shortest | O(E + V log V + K) | Top-K paths (allows repeats) | Yes | ITD detection, efficient top-K |
| Max-flow/min-cut | O(V * E) | Variant boundaries | N/A | INDEL boundary detection |

## 6. Graph Visualization for Debugging

### DOT Format Export

Export the k-mer graph in Graphviz DOT format for visual inspection. This is invaluable for debugging false negatives (why was a variant missed?) and understanding complex graph topologies.

```rust
pub fn to_dot(graph: &KmerGraph) -> String {
    let mut dot = String::from("digraph KmerGraph {\n");
    dot.push_str("  rankdir=LR;\n");
    dot.push_str("  node [shape=box, fontsize=8];\n");

    for (idx, node) in graph.nodes.iter().enumerate() {
        let color = if node.is_reference { "green" }
                    else if node.count > 100 { "orange" }
                    else { "red" };
        let label = format!("{}\\n{}", &node.kmer[..6], node.count);
        dot.push_str(&format!(
            "  n{} [label=\"{}\", style=filled, fillcolor={}];\n",
            idx, label, color
        ));
    }

    for (from, edges) in &graph.edges {
        for edge in edges {
            let style = if edge.is_reference { "solid" } else { "dashed" };
            let penwidth = (1.0 / edge.weight).min(5.0);
            dot.push_str(&format!(
                "  n{} -> n{} [style={}, penwidth={:.1}];\n",
                from, edge.to, style, penwidth
            ));
        }
    }

    dot.push_str("}\n");
    dot
}
```

### Visual Encoding

| Element | Visual Property | Encoding |
|---------|----------------|----------|
| Reference node | Fill color | Green |
| Non-reference node (high count) | Fill color | Orange |
| Non-reference node (low count) | Fill color | Red |
| Node label | Text | First 6 bases + count |
| Reference edge | Line style | Solid |
| Non-reference edge | Line style | Dashed |
| Edge weight | Line width | Inversely proportional to weight |
| Path membership | Node border | Bold if part of a reported variant path |

### Integration with kmerdet CLI

Add a `--debug-graph` flag to the `detect` subcommand:

```
kmerdet detect --jf sample.jf --targets panel/ --debug-graph output/graphs/
```

This writes one DOT file per target (e.g., `NPM1_4ins_exons_10-11utr.dot`). The user can then render with Graphviz:

```bash
dot -Tsvg NPM1_4ins_exons_10-11utr.dot -o NPM1_graph.svg
```

For large graphs (>1000 nodes), the full DOT export is unwieldy. In that case, export only the subgraph within `k` positions of each detected variant, which is typically 50-100 nodes and renders cleanly.

### JSON Graph Format

As an alternative to DOT, export as JSON for consumption by web-based visualization tools (e.g., Cytoscape.js, D3.js):

```json
{
  "nodes": [
    {"id": 0, "kmer": "ACGTACGTACGT...", "count": 5000, "is_reference": true},
    {"id": 1, "kmer": "CGTACGTACGT...", "count": 4950, "is_reference": true}
  ],
  "edges": [
    {"from": 0, "to": 1, "weight": 0.01, "is_reference": true}
  ],
  "paths": [
    {"name": "reference", "nodes": [0, 1, 2, ...]},
    {"name": "variant_SNV_A>T", "nodes": [0, 1, 5, 6, ...]}
  ]
}
```

The JSON format supports the future real-time dashboard feature (see `docs/research/initial/10-future-improvements.md`), where variant graphs could be rendered in a browser.

## References

- Ivanov et al. (2020) -- AStarix: Fast and Optimal Sequence-to-Graph Alignment (A* for genome graphs)
- Yen (1971) -- Finding the K Shortest Loopless Paths in a Network
- Eppstein (1998) -- Finding the k Shortest Paths
- Zerbino & Birney (2008) -- Velvet: de novo short read assembly (Tour Bus bubble collapsing)
- Onodera et al. (2022) -- BubbleGun: enumerating bubbles and superbubbles in genome graphs
- Brankovic et al. (2016) -- Detecting superbubbles in assembly graphs
- Audoux et al. (2017) -- km: original algorithm with Dijkstra-based pathfinding
- Ford & Fulkerson (1956) -- Maximum flow / minimum cut theorem (for flow-based variant boundary detection)
