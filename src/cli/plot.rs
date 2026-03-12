use std::collections::HashMap;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

use crate::variant::{VariantCall, VariantType};

#[derive(clap::Args, Debug)]
pub struct PlotArgs {
    /// Input result file to plot
    #[arg(short, long)]
    pub input: PathBuf,

    /// Chart type to generate
    #[arg(long, default_value = "vaf-histogram")]
    pub chart: ChartType,

    /// Output image file (SVG or PNG)
    #[arg(short, long)]
    pub output: Option<PathBuf>,
}

#[derive(clap::ValueEnum, Debug, Clone, Copy)]
pub enum ChartType {
    /// VAF distribution histogram
    VafHistogram,
    /// Detection rate bar chart
    DetectionBar,
    /// Variant type distribution
    TypeDistribution,
    /// Summary pie chart
    SummaryPie,
}

impl std::fmt::Display for ChartType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::VafHistogram => write!(f, "vaf_histogram"),
            Self::DetectionBar => write!(f, "detection_bar"),
            Self::TypeDistribution => write!(f, "type_distribution"),
            Self::SummaryPie => write!(f, "summary_pie"),
        }
    }
}

/// Parse detection results from a TSV file into VariantCall structs.
fn parse_results(path: &Path) -> Result<Vec<VariantCall>> {
    crate::output::parse_detection_tsv(path)
}

/// Color palette for chart elements.
const COLORS: &[&str] = &[
    "#4a90d9", "#e74c3c", "#2ecc71", "#f39c12", "#9b59b6", "#1abc9c",
    "#e67e22", "#3498db", "#e91e63", "#00bcd4",
];

/// Escape a string for safe use inside SVG text elements.
fn svg_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

// ── SVG generators ──────────────────────────────────────────────────────

/// Generate a VAF histogram SVG for non-reference variant calls.
fn generate_vaf_histogram(calls: &[VariantCall]) -> Result<String> {
    let vafs: Vec<f64> = calls
        .iter()
        .filter(|c| c.variant_type != VariantType::Reference)
        .map(|c| c.rvaf)
        .collect();

    if vafs.is_empty() {
        anyhow::bail!("no non-reference variants to plot");
    }

    // Bin into 20 bins from 0.0 to max_vaf.
    let bin_count = 20usize;
    let max_vaf = vafs.iter().cloned().fold(0.0_f64, f64::max).max(0.01);
    let bin_width = max_vaf / bin_count as f64;
    let mut bins = vec![0u32; bin_count];
    for &vaf in &vafs {
        let bin = ((vaf / bin_width) as usize).min(bin_count - 1);
        bins[bin] += 1;
    }
    let max_count = *bins.iter().max().unwrap_or(&1).max(&1);

    let width = 640;
    let height = 420;
    let margin_left = 65;
    let margin_right = 20;
    let margin_top = 45;
    let margin_bottom = 55;
    let plot_w = width - margin_left - margin_right;
    let plot_h = height - margin_top - margin_bottom;

    let mut svg = format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" style="background: white;">
  <text x="{}" y="28" text-anchor="middle" font-size="16" font-weight="bold" font-family="Arial">VAF Distribution (n={})</text>
"#,
        width,
        height,
        width / 2,
        vafs.len()
    );

    // Draw bars.
    let bar_w = plot_w as f64 / bin_count as f64;
    for (i, &count) in bins.iter().enumerate() {
        let bar_h = count as f64 / max_count as f64 * plot_h as f64;
        let x = margin_left as f64 + i as f64 * bar_w;
        let y = (margin_top + plot_h) as f64 - bar_h;
        svg.push_str(&format!(
            r#"  <rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="{}" stroke="white" stroke-width="1"/>
"#,
            x,
            y,
            bar_w - 1.0,
            bar_h,
            COLORS[0]
        ));
    }

    // X axis line.
    svg.push_str(&format!(
        "  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"#333\" stroke-width=\"1\"/>\n",
        margin_left,
        margin_top + plot_h,
        margin_left + plot_w,
        margin_top + plot_h
    ));

    // Y axis line.
    svg.push_str(&format!(
        "  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"#333\" stroke-width=\"1\"/>\n",
        margin_left, margin_top, margin_left, margin_top + plot_h
    ));

    // X axis labels (5 evenly spaced).
    for i in 0..=4 {
        let val = max_vaf * i as f64 / 4.0;
        let x = margin_left as f64 + i as f64 * plot_w as f64 / 4.0;
        svg.push_str(&format!(
            r#"  <text x="{:.1}" y="{}" text-anchor="middle" font-size="11" font-family="Arial">{:.3}</text>
"#,
            x,
            margin_top + plot_h + 18,
            val
        ));
    }

    // Y axis labels (5 evenly spaced).
    for i in 0..=4 {
        let val = max_count * i / 4;
        let y = (margin_top + plot_h) as f64 - i as f64 * plot_h as f64 / 4.0;
        svg.push_str(&format!(
            r#"  <text x="{}" y="{:.1}" text-anchor="end" font-size="11" font-family="Arial">{}</text>
"#,
            margin_left - 8,
            y + 4.0,
            val
        ));
    }

    // Axis labels.
    svg.push_str(&format!(
        r#"  <text x="{}" y="{}" text-anchor="middle" font-size="13" font-family="Arial">rVAF</text>
"#,
        width / 2,
        height - 5
    ));
    svg.push_str(&format!(
        r#"  <text x="15" y="{}" text-anchor="middle" font-size="13" font-family="Arial" transform="rotate(-90, 15, {})">Count</text>
"#,
        margin_top + plot_h / 2,
        margin_top + plot_h / 2
    ));

    svg.push_str("</svg>\n");
    Ok(svg)
}

/// Generate a detection bar chart SVG showing detected vs reference per target.
fn generate_detection_bar(calls: &[VariantCall]) -> Result<String> {
    if calls.is_empty() {
        anyhow::bail!("no variant calls to plot");
    }

    // Count detected (non-reference) and reference calls per target.
    let mut target_detected: HashMap<String, u32> = HashMap::new();
    let mut target_reference: HashMap<String, u32> = HashMap::new();

    for call in calls {
        if call.variant_type == VariantType::Reference {
            *target_reference.entry(call.target.clone()).or_insert(0) += 1;
        } else {
            *target_detected.entry(call.target.clone()).or_insert(0) += 1;
        }
    }

    // Collect all unique targets in sorted order.
    let mut all_targets: Vec<String> = target_detected
        .keys()
        .chain(target_reference.keys())
        .cloned()
        .collect::<std::collections::HashSet<_>>()
        .into_iter()
        .collect();
    all_targets.sort();

    // Limit to 30 targets for readability; if more, show top by detection count.
    if all_targets.len() > 30 {
        all_targets.sort_by(|a, b| {
            let da = target_detected.get(b).unwrap_or(&0);
            let db_val = target_detected.get(a).unwrap_or(&0);
            da.cmp(db_val)
        });
        all_targets.truncate(30);
        all_targets.sort();
    }

    let n = all_targets.len();
    let max_count = all_targets
        .iter()
        .map(|t| {
            target_detected.get(t).unwrap_or(&0) + target_reference.get(t).unwrap_or(&0)
        })
        .max()
        .unwrap_or(1)
        .max(1);

    // Compute label width: estimate based on longest target name.
    let max_label_len = all_targets.iter().map(|t| t.len()).max().unwrap_or(5);
    let label_space = (max_label_len * 7).min(200).max(80);

    let width = label_space + 500 + 120; // label + plot + legend
    let height = 50 + n * 25 + 50; // top margin + bars + bottom margin
    let margin_left = label_space;
    let margin_top = 45;
    let plot_w = 460;
    let bar_h = 18.0;
    let bar_spacing = 25.0;

    let mut svg = format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" style="background: white;">
  <text x="{}" y="28" text-anchor="middle" font-size="16" font-weight="bold" font-family="Arial">Detection by Target</text>
"#,
        width,
        height,
        width / 2
    );

    for (i, target) in all_targets.iter().enumerate() {
        let detected = *target_detected.get(target).unwrap_or(&0);
        let reference = *target_reference.get(target).unwrap_or(&0);
        let y_base = margin_top as f64 + i as f64 * bar_spacing;

        // Target label.
        let label = if target.len() > 25 {
            format!("{}...", &target[..22])
        } else {
            target.clone()
        };
        svg.push_str(&format!(
            r#"  <text x="{}" y="{:.1}" text-anchor="end" font-size="10" font-family="Arial">{}</text>
"#,
            margin_left - 5,
            y_base + bar_h / 2.0 + 3.5,
            svg_escape(&label)
        ));

        // Detected bar (blue).
        let det_w = detected as f64 / max_count as f64 * plot_w as f64;
        svg.push_str(&format!(
            r#"  <rect x="{}" y="{:.1}" width="{:.1}" height="{:.1}" fill="{}" rx="2"/>
"#,
            margin_left, y_base, det_w, bar_h, COLORS[0]
        ));

        // Reference bar (stacked, gray).
        let ref_w = reference as f64 / max_count as f64 * plot_w as f64;
        svg.push_str(&format!(
            "  <rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" fill=\"#bdc3c7\" rx=\"2\"/>\n",
            margin_left as f64 + det_w,
            y_base,
            ref_w,
            bar_h
        ));

        // Count label.
        let total = detected + reference;
        svg.push_str(&format!(
            "  <text x=\"{:.1}\" y=\"{:.1}\" font-size=\"10\" font-family=\"Arial\" fill=\"#333\">{}</text>\n",
            margin_left as f64 + det_w + ref_w + 4.0,
            y_base + bar_h / 2.0 + 3.5,
            total
        ));
    }

    // Legend.
    let legend_x = margin_left + plot_w + 20;
    let legend_y = margin_top;
    svg.push_str(&format!(
        "  <rect x=\"{}\" y=\"{}\" width=\"12\" height=\"12\" fill=\"{}\"/>\n\
         <text x=\"{}\" y=\"{}\" font-size=\"11\" font-family=\"Arial\">Detected</text>\n\
         <rect x=\"{}\" y=\"{}\" width=\"12\" height=\"12\" fill=\"#bdc3c7\"/>\n\
         <text x=\"{}\" y=\"{}\" font-size=\"11\" font-family=\"Arial\">Reference</text>\n",
        legend_x,
        legend_y,
        COLORS[0],
        legend_x + 16,
        legend_y + 11,
        legend_x,
        legend_y + 22,
        legend_x + 16,
        legend_y + 33
    ));

    svg.push_str("</svg>\n");
    Ok(svg)
}

/// Generate a variant type distribution bar chart SVG.
fn generate_type_distribution(calls: &[VariantCall]) -> Result<String> {
    if calls.is_empty() {
        anyhow::bail!("no variant calls to plot");
    }

    // Count calls per variant type (exclude Reference).
    let mut type_counts: HashMap<VariantType, u32> = HashMap::new();
    let mut ref_count: u32 = 0;
    for call in calls {
        if call.variant_type == VariantType::Reference {
            ref_count += 1;
        } else {
            *type_counts.entry(call.variant_type).or_insert(0) += 1;
        }
    }

    // Ordered list of variant types for consistent display.
    let type_order = [
        VariantType::Substitution,
        VariantType::Insertion,
        VariantType::Deletion,
        VariantType::Itd,
        VariantType::Indel,
    ];

    // Collect types that have non-zero counts.
    let mut bars: Vec<(String, u32, &str)> = Vec::new();
    for (i, vt) in type_order.iter().enumerate() {
        let count = type_counts.get(vt).copied().unwrap_or(0);
        if count > 0 {
            bars.push((vt.to_string(), count, COLORS[i % COLORS.len()]));
        }
    }
    // Also include Reference if present.
    if ref_count > 0 {
        bars.push(("Reference".to_string(), ref_count, "#bdc3c7"));
    }

    if bars.is_empty() {
        anyhow::bail!("no variants to plot");
    }

    let n = bars.len();
    let max_count = bars.iter().map(|(_, c, _)| *c).max().unwrap_or(1).max(1);

    let width = 600;
    let height = 400;
    let margin_left = 80;
    let margin_right = 20;
    let margin_top = 50;
    let margin_bottom = 80;
    let plot_w = width - margin_left - margin_right;
    let plot_h = height - margin_top - margin_bottom;

    let bar_w = (plot_w as f64 / n as f64).min(80.0);
    let total_bar_w = bar_w * n as f64;
    let x_offset = margin_left as f64 + (plot_w as f64 - total_bar_w) / 2.0;

    let mut svg = format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" style="background: white;">
  <text x="{}" y="30" text-anchor="middle" font-size="16" font-weight="bold" font-family="Arial">Variant Type Distribution (n={})</text>
"#,
        width,
        height,
        width / 2,
        calls.len()
    );

    // Axes.
    svg.push_str(&format!(
        "  <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"#333\" stroke-width=\"1\"/>\n\
         <line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"#333\" stroke-width=\"1\"/>\n",
        margin_left,
        margin_top + plot_h,
        margin_left + plot_w,
        margin_top + plot_h,
        margin_left,
        margin_top,
        margin_left,
        margin_top + plot_h
    ));

    // Y axis labels.
    for i in 0..=4 {
        let val = max_count * i / 4;
        let y = (margin_top + plot_h) as f64 - i as f64 * plot_h as f64 / 4.0;
        svg.push_str(&format!(
            r#"  <text x="{}" y="{:.1}" text-anchor="end" font-size="11" font-family="Arial">{}</text>
"#,
            margin_left - 8,
            y + 4.0,
            val
        ));
    }

    // Draw bars and labels.
    for (i, (label, count, color)) in bars.iter().enumerate() {
        let bar_h = *count as f64 / max_count as f64 * plot_h as f64;
        let x = x_offset + i as f64 * bar_w;
        let y = (margin_top + plot_h) as f64 - bar_h;

        svg.push_str(&format!(
            r#"  <rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" fill="{}" stroke="white" stroke-width="1" rx="2"/>
"#,
            x + 2.0,
            y,
            bar_w - 4.0,
            bar_h,
            color
        ));

        // Count label above bar.
        svg.push_str(&format!(
            r#"  <text x="{:.1}" y="{:.1}" text-anchor="middle" font-size="11" font-weight="bold" font-family="Arial">{}</text>
"#,
            x + bar_w / 2.0,
            y - 5.0,
            count
        ));

        // Type label below axis (rotated).
        let label_x = x + bar_w / 2.0;
        let label_y = (margin_top + plot_h + 15) as f64;
        svg.push_str(&format!(
            r#"  <text x="{:.1}" y="{:.1}" text-anchor="start" font-size="11" font-family="Arial" transform="rotate(35, {:.1}, {:.1})">{}</text>
"#,
            label_x, label_y, label_x, label_y, svg_escape(label)
        ));
    }

    // Y axis label.
    svg.push_str(&format!(
        r#"  <text x="15" y="{}" text-anchor="middle" font-size="13" font-family="Arial" transform="rotate(-90, 15, {})">Count</text>
"#,
        margin_top + plot_h / 2,
        margin_top + plot_h / 2
    ));

    svg.push_str("</svg>\n");
    Ok(svg)
}

/// Generate a summary pie chart SVG showing detected vs reference calls.
fn generate_summary_pie(calls: &[VariantCall]) -> Result<String> {
    if calls.is_empty() {
        anyhow::bail!("no variant calls to plot");
    }

    let total = calls.len() as f64;
    let detected = calls
        .iter()
        .filter(|c| c.variant_type != VariantType::Reference)
        .count();
    let reference = calls.len() - detected;

    let width = 500;
    let height = 400;
    let cx = 220.0_f64;
    let cy = 210.0_f64;
    let r = 140.0_f64;

    let mut svg = format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{}" height="{}" style="background: white;">
  <text x="{}" y="30" text-anchor="middle" font-size="16" font-weight="bold" font-family="Arial">Detection Summary (n={})</text>
"#,
        width,
        height,
        width / 2,
        calls.len()
    );

    // Build slices: (label, count, color).
    let slices: Vec<(&str, usize, &str)> = vec![
        ("Detected", detected, COLORS[0]),
        ("Reference", reference, "#bdc3c7"),
    ];

    let mut start_angle: f64 = -std::f64::consts::FRAC_PI_2; // Start from top (12 o'clock).
    for (_label, count, color) in &slices {
        if *count == 0 {
            continue;
        }
        let fraction = *count as f64 / total;
        let sweep_angle = fraction * 2.0 * std::f64::consts::PI;
        let end_angle = start_angle + sweep_angle;

        let large_arc = if sweep_angle > std::f64::consts::PI {
            1
        } else {
            0
        };

        let x1 = cx + r * start_angle.cos();
        let y1 = cy + r * start_angle.sin();
        let x2 = cx + r * end_angle.cos();
        let y2 = cy + r * end_angle.sin();

        // Handle full circle (single slice = 100%).
        if fraction >= 0.999 {
            svg.push_str(&format!(
                r#"  <circle cx="{:.1}" cy="{:.1}" r="{:.1}" fill="{}"/>
"#,
                cx, cy, r, color
            ));
        } else {
            svg.push_str(&format!(
                r#"  <path d="M {:.1} {:.1} L {:.2} {:.2} A {:.1} {:.1} 0 {} 1 {:.2} {:.2} Z" fill="{}" stroke="white" stroke-width="2"/>
"#,
                cx, cy, x1, y1, r, r, large_arc, x2, y2, color
            ));
        }

        // Label at the midpoint of the arc.
        let mid_angle = start_angle + sweep_angle / 2.0;
        let label_r = r * 0.65;
        let lx = cx + label_r * mid_angle.cos();
        let ly = cy + label_r * mid_angle.sin();
        let pct = fraction * 100.0;
        if pct >= 3.0 {
            svg.push_str(&format!(
                r#"  <text x="{:.1}" y="{:.1}" text-anchor="middle" font-size="12" font-weight="bold" font-family="Arial" fill="white">{:.1}%</text>
"#,
                lx, ly, pct
            ));
        }

        start_angle = end_angle;
    }

    // Legend on the right side.
    let legend_x = 390;
    let legend_y = 80;
    for (i, (label, count, color)) in slices.iter().enumerate() {
        let ly = legend_y + i * 28;
        svg.push_str(&format!(
            r#"  <rect x="{}" y="{}" width="14" height="14" fill="{}" rx="2"/>
  <text x="{}" y="{}" font-size="12" font-family="Arial">{} ({})</text>
"#,
            legend_x,
            ly,
            color,
            legend_x + 20,
            ly + 12,
            label,
            count
        ));
    }

    svg.push_str("</svg>\n");
    Ok(svg)
}

pub fn run(args: PlotArgs, _global: &super::GlobalOptions) -> Result<()> {
    // 1. Read input file (TSV detection results).
    let calls = parse_results(&args.input)?;
    tracing::info!("Loaded {} variant calls from {}", calls.len(), args.input.display());

    if calls.is_empty() {
        anyhow::bail!("no variant calls found in {}", args.input.display());
    }

    // 2. Generate chart based on type.
    let svg = match args.chart {
        ChartType::VafHistogram => generate_vaf_histogram(&calls)?,
        ChartType::DetectionBar => generate_detection_bar(&calls)?,
        ChartType::TypeDistribution => generate_type_distribution(&calls)?,
        ChartType::SummaryPie => generate_summary_pie(&calls)?,
    };

    // 3. Determine output path.
    let output_path = args.output.unwrap_or_else(|| {
        PathBuf::from(format!("{}.svg", args.chart))
    });

    // 4. Write SVG output.
    std::fs::write(&output_path, &svg)
        .with_context(|| format!("writing SVG to {}", output_path.display()))?;
    tracing::info!("Plot written to {} ({} bytes)", output_path.display(), svg.len());

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_calls() -> Vec<VariantCall> {
        vec![
            VariantCall {
                sample: "s1".to_string(),
                target: "targetA".to_string(),
                variant_type: VariantType::Reference,
                variant_name: String::new(),
                rvaf: 1.0,
                expression: 100.0,
                min_coverage: 50,
                path_score: 50,
                start_kmer_count: 50,
                ref_sequence: "ACGT".to_string(),
                alt_sequence: "ACGT".to_string(),
                info: "vs_ref".to_string(),
                chrom: None,
                pos: None,
                ref_allele: None,
                alt_allele: None,
                pvalue: None,
                qual: None,
                ci_lower: None,
                ci_upper: None,
            },
            VariantCall {
                sample: "s1".to_string(),
                target: "targetA".to_string(),
                variant_type: VariantType::Substitution,
                variant_name: "4:A/T:4".to_string(),
                rvaf: 0.15,
                expression: 15.0,
                min_coverage: 10,
                path_score: 10,
                start_kmer_count: 10,
                ref_sequence: "ACGT".to_string(),
                alt_sequence: "TCGT".to_string(),
                info: "vs_ref".to_string(),
                chrom: Some("chr1".to_string()),
                pos: Some(12345),
                ref_allele: Some("A".to_string()),
                alt_allele: Some("T".to_string()),
                pvalue: None,
                qual: None,
                ci_lower: None,
                ci_upper: None,
            },
            VariantCall {
                sample: "s1".to_string(),
                target: "targetB".to_string(),
                variant_type: VariantType::Insertion,
                variant_name: "10:A/AGGG:12".to_string(),
                rvaf: 0.08,
                expression: 8.0,
                min_coverage: 5,
                path_score: 5,
                start_kmer_count: 5,
                ref_sequence: "ACGT".to_string(),
                alt_sequence: "AGGGCGT".to_string(),
                info: "vs_ref".to_string(),
                chrom: Some("chr2".to_string()),
                pos: Some(67890),
                ref_allele: Some("A".to_string()),
                alt_allele: Some("AGGG".to_string()),
                pvalue: None,
                qual: None,
                ci_lower: None,
                ci_upper: None,
            },
            VariantCall {
                sample: "s1".to_string(),
                target: "targetB".to_string(),
                variant_type: VariantType::Reference,
                variant_name: String::new(),
                rvaf: 1.0,
                expression: 80.0,
                min_coverage: 40,
                path_score: 40,
                start_kmer_count: 40,
                ref_sequence: "ACGT".to_string(),
                alt_sequence: "ACGT".to_string(),
                info: "vs_ref".to_string(),
                chrom: None,
                pos: None,
                ref_allele: None,
                alt_allele: None,
                pvalue: None,
                qual: None,
                ci_lower: None,
                ci_upper: None,
            },
            VariantCall {
                sample: "s1".to_string(),
                target: "targetC".to_string(),
                variant_type: VariantType::Deletion,
                variant_name: "5:GT/:7".to_string(),
                rvaf: 0.22,
                expression: 22.0,
                min_coverage: 12,
                path_score: 12,
                start_kmer_count: 12,
                ref_sequence: "ACGTACGT".to_string(),
                alt_sequence: "ACACGT".to_string(),
                info: "vs_ref".to_string(),
                chrom: Some("chr3".to_string()),
                pos: Some(11111),
                ref_allele: Some("GT".to_string()),
                alt_allele: Some("".to_string()),
                pvalue: None,
                qual: None,
                ci_lower: None,
                ci_upper: None,
            },
        ]
    }

    #[test]
    fn test_vaf_histogram_svg() {
        let calls = make_test_calls();
        let svg = generate_vaf_histogram(&calls).unwrap();
        assert!(svg.starts_with("<svg"));
        assert!(svg.contains("VAF Distribution"));
        assert!(svg.contains("</svg>"));
    }

    #[test]
    fn test_detection_bar_svg() {
        let calls = make_test_calls();
        let svg = generate_detection_bar(&calls).unwrap();
        assert!(svg.starts_with("<svg"));
        assert!(svg.contains("Detection by Target"));
        assert!(svg.contains("</svg>"));
    }

    #[test]
    fn test_type_distribution_svg() {
        let calls = make_test_calls();
        let svg = generate_type_distribution(&calls).unwrap();
        assert!(svg.starts_with("<svg"));
        assert!(svg.contains("Variant Type Distribution"));
        assert!(svg.contains("</svg>"));
    }

    #[test]
    fn test_summary_pie_svg() {
        let calls = make_test_calls();
        let svg = generate_summary_pie(&calls).unwrap();
        assert!(svg.starts_with("<svg"));
        assert!(svg.contains("Detection Summary"));
        assert!(svg.contains("</svg>"));
    }

    #[test]
    fn test_chart_type_display() {
        assert_eq!(ChartType::VafHistogram.to_string(), "vaf_histogram");
        assert_eq!(ChartType::DetectionBar.to_string(), "detection_bar");
        assert_eq!(ChartType::TypeDistribution.to_string(), "type_distribution");
        assert_eq!(ChartType::SummaryPie.to_string(), "summary_pie");
    }

    #[test]
    fn test_vaf_histogram_empty_non_reference() {
        // Only reference calls -- should return an error.
        let calls = vec![VariantCall {
            sample: "s1".to_string(),
            target: "t1".to_string(),
            variant_type: VariantType::Reference,
            variant_name: String::new(),
            rvaf: 1.0,
            expression: 100.0,
            min_coverage: 50,
            path_score: 50,
            start_kmer_count: 50,
            ref_sequence: "ACGT".to_string(),
            alt_sequence: "ACGT".to_string(),
            info: String::new(),
            chrom: None,
            pos: None,
            ref_allele: None,
            alt_allele: None,
            pvalue: None,
            qual: None,
            ci_lower: None,
            ci_upper: None,
        }];
        let result = generate_vaf_histogram(&calls);
        assert!(result.is_err());
    }

    #[test]
    fn test_svg_escape() {
        assert_eq!(svg_escape("A<B>C&D\"E"), "A&lt;B&gt;C&amp;D&quot;E");
        assert_eq!(svg_escape("normal"), "normal");
    }
}
