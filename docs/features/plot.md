# Feature: `plot` -- Visualization

## What It Does

The `plot` subcommand generates publication-quality visualizations from variant detection
results. It produces standard plot types used in clinical reporting and research analysis:
VAF histograms, detection summary charts, variant type distributions, and summary pie
charts. An additional longitudinal VAF plot tracks variant allele frequencies across
serial timepoints for MRD monitoring.

```
kmerdet plot -i filtered.tsv -o plots/
kmerdet plot -i filtered.tsv --type vaf-histogram --format png -o vaf_dist.png
kmerdet plot --longitudinal -i T1.tsv T2.tsv T3.tsv -o vaf_trend.svg \
             --timepoint-labels "Baseline,Week4,Week8"
```

### Input Specification

| Input | Flag | Format | Required |
|-------|------|--------|----------|
| Result file(s) | `-i` / `--input` | TSV/CSV from detect or filter | Yes |
| Output path | `-o` / `--output` | File path or directory | Yes |
| Plot type | `--type` | Specific plot or "all" | No (all) |
| Format | `--format` | svg / png | No (svg) |

### Output Specification

| Plot | Default filename | Description |
|------|-----------------|-------------|
| VAF histogram | `vaf_histogram.svg` | Distribution of rVAF values |
| Detection summary | `detection_summary.svg` | Detected vs. not-detected per target |
| Type distribution | `type_distribution.svg` | Variant type breakdown |
| Summary pie | `summary_pie.svg` | Overall detection rate |
| Longitudinal VAF | `vaf_longitudinal.svg` | VAF trend across timepoints |

---

## Why It Matters

### Clinical Reporting

Visual summaries are a critical component of clinical reports for liquid biopsy results.
Oncologists reviewing MRD monitoring data need to see trends at a glance: Is the patient
responding to treatment? Are any tracked mutations increasing in VAF? A VAF trend plot
showing declining values over time communicates treatment response far more effectively
than a table of numbers.

The thesis pipeline relied on matplotlib via kmtools plot (available only in the
development branch). Moving visualization into the core Rust binary eliminates the
Python/matplotlib dependency and ensures that plots are available in every deployment,
including minimal clinical environments without Python.

### Research and Validation

During pipeline development and validation, plots reveal patterns invisible in tabular
data. The VAF histogram exposes bimodal distributions (separating true variants from
noise), the type distribution reveals whether INDEL sensitivity improvements are
working, and the detection summary identifies targets that consistently fail.

---

## Plot Types

### VAF Histogram

A histogram of rVAF values across all detected variants (excluding Reference calls).
This is the most commonly used visualization for ctDNA detection results.

**Configuration**:

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Bin width | `--bin-width` | 0.005 | Width of each histogram bin |
| Min VAF | `--vaf-min` | 0.0 | Left axis limit |
| Max VAF | `--vaf-max` | auto | Right axis limit (auto-scaled to data) |
| Log scale | `--log-y` | false | Logarithmic y-axis |
| Color | `--color` | "#4C72B0" | Bar fill color (CSS color string) |

**Features**:
- Vertical dashed line at median VAF
- Annotation showing mean, median, and count
- Optional overlay of a second dataset (for before/after comparison)
- X-axis labeled "Relative VAF" with appropriate decimal formatting
- Y-axis labeled "Number of variants"

**Clinical interpretation**: A tight cluster of VAF values near 0.01-0.05 suggests
stable low-level ctDNA (MRD positive). A broad distribution spanning 0.01-0.50 suggests
heterogeneous tumor shedding. No variants above threshold suggests MRD negative.

### Detection Summary Bar Chart

A grouped bar chart showing, for each target, whether a variant was detected. Targets
are listed on the y-axis (sorted alphabetically or by detection status), with bars
indicating detected (colored) vs. not detected (gray).

**Configuration**:

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Sort order | `--sort` | name | name, status, or vaf |
| Show VAF | `--show-vaf` | false | Color bars by rVAF value |
| Max targets | `--max-targets` | 100 | Truncate display for large panels |

**Features**:
- Color coding: green for detected, gray for not detected, orange for failed
- Optional rVAF annotation next to each bar
- Summary counts at the top: "38/50 detected (76%)"

### Variant Type Distribution

A bar chart showing the count of each variant type:

| Type | Color |
|------|-------|
| Substitution | Blue |
| Insertion | Green |
| Deletion | Red |
| ITD | Purple |
| Complex | Orange |
| Reference | Gray |

**Features**:
- Counts labeled on each bar
- Optional percentage labels
- Optional grouping by target (stacked bar chart)

### Summary Pie Chart

A pie chart showing the overall distribution of variant call outcomes:

| Slice | Description |
|-------|-------------|
| Detected (PASS) | Variants passing all filters |
| Filtered | Variants detected but removed by filtering |
| Not detected | Targets with Reference-only calls |
| Failed | Targets that produced no result |

**Features**:
- Percentage labels on each slice
- Total count in the center
- Color scheme consistent with the detection summary chart

### Longitudinal VAF Plot

A line plot tracking VAF values across serial timepoints. This is the most clinically
important visualization for MRD monitoring.

```
kmerdet plot --longitudinal -i T1.tsv T2.tsv T3.tsv T4.tsv \
             --timepoint-labels "Baseline,Week4,Week8,Week12" \
             -o vaf_trend.svg
```

**Configuration**:

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Timepoint labels | `--timepoint-labels` | T1,T2,...,Tn | X-axis labels |
| Variants to track | `--track` | all | Specific variant names or "all" |
| Detection threshold | `--threshold` | 0.001 | Horizontal line at detection limit |
| Log scale | `--log-y` | true | Logarithmic y-axis (typical for VAF) |
| Show CI | `--show-ci` | false | Display confidence intervals as ribbons |

**Features**:
- One line per tracked variant, colored distinctly
- Points where a variant drops below detection (VAF = 0) shown as open circles at
  the detection threshold
- Horizontal dashed line at the detection threshold (0.1% VAF)
- Legend mapping colors to variant names
- Optional trend line (linear regression on log-VAF)
- Optional confidence interval ribbons when bootstrap CIs are available

**Clinical interpretation**: Declining VAF lines = treatment response. Rising VAF lines
after initial decline = molecular relapse (requires clinical action). Stable low VAF =
stable MRD (continue monitoring).

---

## Output Formats

### SVG (Default)

Scalable Vector Graphics. Advantages:
- Lossless at any zoom level
- Editable in vector graphics editors (Inkscape, Illustrator)
- Small file size for simple charts
- Embeddable in HTML reports

### PNG

Rasterized output for inclusion in presentations and static reports. Rendered at
300 DPI by default (configurable via `--dpi`).

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| DPI | `--dpi` | 300 | Dots per inch for PNG output |
| Width | `--width` | 800 | Image width in pixels |
| Height | `--height` | 600 | Image height in pixels |

---

## Customization

### Global Options

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Title | `--title` | Auto-generated | Plot title |
| Font size | `--font-size` | 14 | Base font size in points |
| Color scheme | `--colors` | default | Named palette or custom hex list |
| Background | `--background` | white | Background color |
| Dimensions | `--width`, `--height` | 800x600 | Canvas dimensions |

### Named Color Palettes

| Palette | Description | Source |
|---------|-------------|--------|
| default | Blue-focused, colorblind-safe | Seaborn "muted" inspired |
| clinical | High contrast, print-friendly | Custom for clinical reports |
| viridis | Perceptually uniform | Matplotlib viridis |
| grayscale | Black/white/gray only | For B&W printing |

---

## Implementation

### Plotters Crate

All visualizations are implemented using the `plotters` crate, chosen for:
- Pure Rust (no system dependencies)
- SVG and bitmap backends
- Sufficient API for all required chart types
- Active maintenance and Rust ecosystem integration

The rendering pipeline:
1. Parse input data (reuse the same TSV/CSV parser as `stats`)
2. Compute plot-specific derived data (bins, percentages, aggregations)
3. Build a plotters chart with appropriate axes, series, and annotations
4. Render to the selected backend (SVG file or PNG bitmap)

### Separation of Concerns

Plot logic is separated into data preparation (computing what to plot) and rendering
(how to draw it). Data preparation is shared with the `stats` subcommand where
applicable (e.g., VAF percentiles, type counts). Rendering is specific to the plot
subcommand.

---

## Research Backing

### Clinical Workflow Need

The thesis describes a clinical workflow where rapid k-mer detection produces same-day
preliminary results. Visual summaries are part of this reporting: oncologists reviewing
MRD results expect to see a VAF trend plot alongside the tabular results. Commercial
ctDNA platforms (Guardant360, FoundationOne Liquid CDx) include visualization in their
reports as standard practice.

### Longitudinal Monitoring

The thesis proposed using k-mer detection for longitudinal VAF tracking across serial
blood draws. The longitudinal VAF plot directly implements this: each timepoint's
results are plotted as a point, and the trend reveals treatment response or relapse.
The confidence-metrics.md research document notes that bootstrap CIs on rVAF would
enable plotting uncertainty ribbons, distinguishing true VAF changes from measurement
noise.

### Sensitivity Visualization

The sensitivity-landscape.md document contains comparative sensitivity data across
VAF ranges and variant types. Generating type-specific VAF histograms from kmerdet
results enables direct comparison with the thesis benchmarks and identification of
sensitivity improvements from algorithmic changes.

---

## Design Considerations

### Headless Operation

The plot subcommand runs without a display server (no X11/Wayland required). This is
essential for clinical pipelines running on headless servers and HPC clusters. The
plotters crate renders directly to file without any GUI dependency.

### Deterministic Output

SVG output is deterministic: the same input data produces byte-identical SVG files.
PNG output may vary slightly across platforms due to font rendering differences, but
the visual content is identical.

### Large Panels

For panels with >100 targets, the detection summary chart becomes unwieldy. The
`--max-targets` flag limits the display, with a note indicating truncation. Targets
are prioritized by clinical relevance (detected variants first, then alphabetical).

---

## Acceptance Criteria

### Unit Tests

- [ ] VAF histogram: correct bin assignments for known rVAF values
- [ ] Detection summary: correct detected/not-detected counts
- [ ] Type distribution: correct per-type counts match stats output
- [ ] Summary pie: slice percentages sum to 100%
- [ ] Longitudinal: correct mapping of timepoints to data points

### Integration Tests

- [ ] SVG output is valid XML (parseable by an XML parser)
- [ ] PNG output is valid PNG (readable by image library)
- [ ] All plot types render without errors on test dataset
- [ ] Default dimensions produce readable labels (no text overlap)
- [ ] Custom color schemes are applied correctly
- [ ] Longitudinal plot with missing timepoints (variant not detected) handles gaps

### Visual Verification

- [ ] VAF histogram: bars are proportional to counts, median line is correct
- [ ] Detection summary: all targets visible, colors distinguish detected from not
- [ ] Type distribution: bar heights match expected counts
- [ ] Longitudinal: lines connect correct timepoints, log scale is legible
- [ ] All plots: axis labels, titles, and legends are readable at default size
- [ ] Grayscale palette: all elements distinguishable in B&W print
