//! Dependency-free SVG figures for protstruct.
//!
//! Modes:
//!   protstruct-plot workflow --output workflow.svg
//!   protstruct-plot consistency --input consistency.tsv --output consistency.svg
//!   protstruct-plot mtba --input mtba.tsv --output mtba.svg
//!
//! No external crates are used. SVG text uses Times New Roman + SimSun fallback.

use std::env;
use std::fmt::Write as FmtWrite;
use std::fs;
use std::process;

const FONT: &str = "Times New Roman, SimSun, Songti SC, Noto Serif CJK SC, serif";
const INK: &str = "#334155";
const SLATE: &str = "#5878A6";
const BLUE: &str = "#A7C7E7";
const PINK: &str = "#F4A7AE";
const PALE: &str = "#CFEAF6";
const PEACH: &str = "#FBD9C9";
const MINT: &str = "#A6D0C6";
const GRAY: &str = "#DADDE3";
const BG: &str = "#F7FAFC";

#[derive(Default, Clone)]
struct Row {
    id: String,
    ss: String,
    total: f64,
    reference: bool,
}

struct Consistency {
    reference: String,
    max_total: f64,
    width: usize,
    agreement: Vec<f64>,
    rows: Vec<Row>,
    disulfides: Vec<(usize, usize)>,
    salt: Vec<(usize, usize)>,
}

#[derive(Clone)]
struct Composition {
    label: String,
    h: f64,
    e: f64,
    c: f64,
    total: f64,
    ratio: f64,
}

fn tr<'a>(lang: &str, zh: &'a str, en: &'a str) -> &'a str {
    if lang == "zh" { zh } else { en }
}

fn escape(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

fn text(out: &mut String, x: f64, y: f64, size: f64, anchor: &str, fill: &str, weight: &str, value: &str) {
    let _ = write!(
        out,
        "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"{}\" font-size=\"{:.1}\" text-anchor=\"{}\" fill=\"{}\" font-weight=\"{}\">{}</text>\n",
        x, y, FONT, size, anchor, fill, weight, escape(value)
    );
}

fn text_rot(out: &mut String, x: f64, y: f64, size: f64, anchor: &str, fill: &str, value: &str, angle: f64) {
    let _ = write!(
        out,
        "<text x=\"{:.1}\" y=\"{:.1}\" transform=\"rotate({:.1} {:.1} {:.1})\" font-family=\"{}\" font-size=\"{:.1}\" text-anchor=\"{}\" fill=\"{}\">{}</text>\n",
        x, y, angle, x, y, FONT, size, anchor, fill, escape(value)
    );
}

fn rect(out: &mut String, x: f64, y: f64, w: f64, h: f64, fill: &str, stroke: &str, width: f64, rx: f64) {
    let _ = write!(
        out,
        "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" rx=\"{:.1}\" fill=\"{}\" stroke=\"{}\" stroke-width=\"{:.1}\"/>\n",
        x, y, w, h, rx, fill, stroke, width
    );
}

fn line(out: &mut String, x1: f64, y1: f64, x2: f64, y2: f64, stroke: &str, width: f64) {
    let _ = write!(
        out,
        "<line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"{}\" stroke-width=\"{:.1}\"/>\n",
        x1, y1, x2, y2, stroke, width
    );
}

fn consistency_color(char: char) -> &'static str {
    match char {
        'H' => PINK,
        'E' => BLUE,
        'C' => GRAY,
        _ => "#FFFFFF",
    }
}

fn parse_consistency(path: &str) -> Result<Consistency, String> {
    let content = fs::read_to_string(path).map_err(|error| format!("read {}: {}", path, error))?;
    let mut data = Consistency {
        reference: String::new(),
        max_total: 1.0,
        width: 1,
        agreement: Vec::new(),
        rows: Vec::new(),
        disulfides: Vec::new(),
        salt: Vec::new(),
    };
    for line in content.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        match fields.first().copied().unwrap_or("") {
            "REF" => data.reference = fields.get(1).unwrap_or(&"").to_string(),
            "MAX" => data.max_total = fields.get(1).unwrap_or(&"1").parse().unwrap_or(1.0),
            "WIDTH" => data.width = fields.get(1).unwrap_or(&"1").parse().unwrap_or(1),
            "AGREE" => {
                data.agreement = fields
                    .get(1)
                    .unwrap_or(&"")
                    .split(',')
                    .filter_map(|value| value.parse().ok())
                    .collect();
            }
            "SS" => {
                if fields.len() >= 3 {
                    data.rows.push(Row {
                        id: fields[1].to_string(),
                        ss: fields[2].to_string(),
                        total: 0.0,
                        reference: false,
                    });
                }
            }
            "SCORE" => {
                if fields.len() >= 4 {
                    let id = fields[1];
                    let total = fields[2].parse().unwrap_or(0.0);
                    let reference = fields[3] == "1";
                    if let Some(row) = data.rows.iter_mut().find(|row| row.id == id) {
                        row.total = total;
                        row.reference = reference;
                    }
                }
            }
            "DS" | "SB" => {
                if fields.len() >= 3 {
                    let row = fields[1].parse().unwrap_or(0);
                    let column = fields[2].parse().unwrap_or(0);
                    if fields[0] == "DS" {
                        data.disulfides.push((row, column));
                    } else {
                        data.salt.push((row, column));
                    }
                }
            }
            _ => {}
        }
    }
    if data.rows.is_empty() {
        return Err("consistency data has no SS rows".to_string());
    }
    if data.width == 0 {
        data.width = data.rows.iter().map(|row| row.ss.chars().count()).max().unwrap_or(1);
    }
    Ok(data)
}

fn consistency_svg(data: &Consistency, lang: &str) -> String {
    let rows = data.rows.len().max(1);
    let width = data.width.max(1);
    let height = 360.0 + rows as f64 * 58.0;
    let left = 170.0;
    let map_w = 1080.0;
    let map_y = 210.0;
    let row_h = 58.0;
    let map_h = rows as f64 * row_h;
    let score_x = left + map_w + 35.0;
    let score_w = 215.0;
    let cell = map_w / width as f64;
    let mut svg = String::new();
    let _ = writeln!(
        svg,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1600\" height=\"{:.0}\" viewBox=\"0 0 1600 {:.0}\" shape-rendering=\"crispEdges\">",
        height, height
    );
    rect(&mut svg, 0.0, 0.0, 1600.0, height, "#FFFFFF", "#FFFFFF", 0.0, 0.0);
    let title = if lang == "zh" {
        format!("多序列二级结构一致性图（参考序列：{}）", data.reference)
    } else {
        format!("Multi-sequence secondary-structure consistency (reference: {})", data.reference)
    };
    text(&mut svg, 800.0, 30.0, 24.0, "middle", INK, "bold", &title);
    text(&mut svg, 800.0, 56.0, 13.0, "middle", SLATE, "normal", tr(lang, "逐列比较：二级结构相同 +1，不同 0；二硫键/盐桥对应 +0.5；参考序列为第一条序列", "Column comparison: SS match +1, mismatch 0; feature-pair match +0.5; first sequence is the reference"));

    // Agreement panel.
    let agreement_y = 82.0;
    let agreement_h = 88.0;
    line(&mut svg, left, agreement_y + agreement_h, left + map_w, agreement_y + agreement_h, BLUE, 1.0);
    for (index, value) in data.agreement.iter().take(width).enumerate() {
        let h = value.clamp(0.0, 1.0) * agreement_h;
        rect(&mut svg, left + index as f64 * cell - 0.15, agreement_y + agreement_h - h, cell + 0.3, h, BLUE, BLUE, 0.0, 0.0);
    }
    text(&mut svg, left - 12.0, agreement_y + 8.0, 12.0, "end", INK, "normal", "100%");
    text(&mut svg, left - 12.0, agreement_y + agreement_h * 0.55, 12.0, "end", INK, "normal", "50%");
    text(&mut svg, left - 12.0, agreement_y + agreement_h, 12.0, "end", INK, "normal", "0%");
    text_rot(&mut svg, left - 48.0, agreement_y + agreement_h * 0.55, 12.0, "middle", INK, tr(lang, "逐列一致性", "Column agreement"), -90.0);

    // Heatmap.
    for row_index in 0..rows {
        let y = map_y + row_index as f64 * row_h;
        rect(&mut svg, left, y, map_w, row_h, "#FFFFFF", "#FFFFFF", 0.0, 0.0);
        for (column, residue) in data.rows[row_index].ss.chars().take(width).enumerate() {
            rect(
                &mut svg,
                left + column as f64 * cell - 0.15,
                y,
                cell + 0.3,
                row_h - 1.0,
                consistency_color(residue),
                consistency_color(residue),
                0.0,
                0.0,
            );
        }
        let label = if data.rows[row_index].reference {
            if lang == "en" {
                format!("{} (reference)", data.rows[row_index].id)
            } else {
                format!("{}（参考）", data.rows[row_index].id)
            }
        } else {
            data.rows[row_index].id.clone()
        };
        text(&mut svg, left - 14.0, y + row_h * 0.58, 13.0, "end", INK, "normal", &label);
        line(&mut svg, left, y + row_h, left + map_w, y + row_h, "#FFFFFF", 1.0);
    }

    // Feature markers.
    for (row, column) in &data.disulfides {
        if *row < rows && *column < width {
            let _ = write!(
                svg,
                "<circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"3.3\" fill=\"{}\"/>\n",
                left + (*column as f64 + 0.5) * cell,
                map_y + (*row as f64 + 0.5) * row_h,
                SLATE
            );
        }
    }
    for (row, column) in &data.salt {
        if *row < rows && *column < width {
            let _ = write!(
                svg,
                "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"5.5\" height=\"5.5\" fill=\"{}\" stroke=\"{}\" stroke-width=\"1\"/>\n",
                left + (*column as f64 + 0.5) * cell - 2.7,
                map_y + (*row as f64 + 0.5) * row_h - 2.7,
                PEACH,
                SLATE
            );
        }
    }
    for tick in (0..width).step_by(if width <= 60 { 10 } else if width <= 200 { 20 } else { 50 }) {
        text(&mut svg, left + (tick as f64 + 0.5) * cell, map_y + map_h + 20.0, 11.0, "middle", INK, "normal", &(tick + 1).to_string());
    }
    text(&mut svg, left + map_w / 2.0, map_y + map_h + 42.0, 12.0, "middle", INK, "normal", tr(lang, "比对位置（列）", "Aligned position (column)"));

    // Score panel.
    let score_max = data.max_total.max(1.0);
    line(&mut svg, score_x, map_y + map_h, score_x + score_w, map_y + map_h, SLATE, 1.0);
    line(&mut svg, score_x, map_y, score_x, map_y + map_h, SLATE, 1.0);
    for row_index in 0..rows {
        let bar_w = (data.rows[row_index].total / score_max).clamp(0.0, 1.0) * score_w;
        let color = if data.rows[row_index].reference { SLATE } else { BLUE };
        rect(&mut svg, score_x, map_y + row_index as f64 * row_h + 12.0, bar_w, row_h - 24.0, color, color, 0.0, 0.0);
        text(&mut svg, score_x + bar_w + 7.0, map_y + (row_index as f64 + 0.5) * row_h + 4.0, 11.0, "start", INK, "normal", &format!("{:.1} ({:.1}%)", data.rows[row_index].total, data.rows[row_index].total / score_max * 100.0));
    }
    text(&mut svg, score_x + score_w / 2.0, map_y + map_h + 42.0, 12.0, "middle", INK, "normal", tr(lang, "总分 / 总分一致率", "Total / score ratio"));

    // Legend.
    let legend_y = map_y + map_h + 78.0;
    let legend = [
        (PINK, tr(lang, "α-螺旋 (H)", "α-helix (H)"), "rect"),
        (BLUE, tr(lang, "β-折叠 (E)", "β-sheet (E)"), "rect"),
        (GRAY, tr(lang, "coil/卷曲 (C)", "Coil (C)"), "rect"),
        (SLATE, tr(lang, "二硫键残基", "Disulfide"), "circle"),
        (PEACH, tr(lang, "盐桥残基", "Salt bridge"), "square"),
    ];
    let mut x = 310.0;
    for (color, label, kind) in legend {
        if kind == "circle" {
            let _ = write!(svg, "<circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"4\" fill=\"{}\"/>\n", x, legend_y, color);
        } else if kind == "square" {
            let _ = write!(svg, "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"8\" height=\"8\" fill=\"{}\" stroke=\"{}\" stroke-width=\"1\"/>\n", x - 4.0, legend_y - 4.0, color, SLATE);
        } else {
            rect(&mut svg, x - 6.0, legend_y - 6.0, 12.0, 12.0, color, "#FFFFFF", 0.5, 0.0);
        }
        text(&mut svg, x + 12.0, legend_y + 4.0, 13.0, "start", INK, "normal", label);
        x += 205.0;
    }
    svg.push_str("</svg>\n");
    svg
}

fn parse_mtba(path: &str) -> Result<(Vec<Composition>, Vec<Composition>), String> {
    let content = fs::read_to_string(path).map_err(|error| format!("read {}: {}", path, error))?;
    let mut paralogs = Vec::new();
    let mut orthologs = Vec::new();
    for line in content.lines() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 7 {
            continue;
        }
        let item = Composition {
            label: fields[1].to_string(),
            h: fields[2].parse().unwrap_or(0.0),
            e: fields[3].parse().unwrap_or(0.0),
            c: fields[4].parse().unwrap_or(0.0),
            total: fields[5].parse().unwrap_or(0.0),
            ratio: fields[6].parse().unwrap_or(0.0),
        };
        if fields[0] == "PARA" {
            paralogs.push(item);
        } else if fields[0] == "ORTH" {
            orthologs.push(item);
        }
    }
    Ok((paralogs, orthologs))
}

fn composition_panel(
    out: &mut String,
    x: f64,
    y: f64,
    width: f64,
    height: f64,
    title: &str,
    data: &[Composition],
    lang: &str,
) {
    text(out, x, y - 12.0, 18.0, "start", INK, "bold", title);
    let baseline = y + height;
    let slot = width / data.len().max(1) as f64;
    let bar_w = (slot * 0.58).min(72.0);
    for (index, item) in data.iter().enumerate() {
        let cx = x + slot * (index as f64 + 0.5);
        let h_h = item.h / 100.0 * height;
        let e_h = item.e / 100.0 * height;
        let c_h = item.c / 100.0 * height;
        rect(out, cx - bar_w / 2.0, baseline - h_h, bar_w, h_h, PINK, "#FFFFFF", 0.8, 0.0);
        rect(out, cx - bar_w / 2.0, baseline - h_h - e_h, bar_w, e_h, BLUE, "#FFFFFF", 0.8, 0.0);
        rect(out, cx - bar_w / 2.0, baseline - h_h - e_h - c_h, bar_w, c_h, GRAY, "#FFFFFF", 0.8, 0.0);
        text(out, cx, baseline - h_h / 2.0 + 4.0, 12.0, "middle", INK, "bold", &format!("{:.1}", item.h));
        text_rot(out, cx, baseline + 18.0, 11.0, "end", INK, &item.label, 25.0);
    }
    line(out, x, baseline, x + width, baseline, SLATE, 1.0);
    text_rot(out, x - 35.0, y + height / 2.0, 12.0, "middle", INK, tr(lang, "占比 (%)", "Composition (%)"), -90.0);
    for value in [0, 50, 100] {
        let py = baseline - value as f64 / 100.0 * height;
        line(out, x - 6.0, py, x, py, SLATE, 1.0);
        text(out, x - 10.0, py + 4.0, 11.0, "end", INK, "normal", &value.to_string());
    }
}

fn score_panel(
    out: &mut String,
    x: f64,
    y: f64,
    width: f64,
    height: f64,
    title: &str,
    data: &[Composition],
    reference: &str,
    lang: &str,
) {
    text(out, x, y - 12.0, 18.0, "start", INK, "bold", title);
    let baseline = y + height;
    let slot = width / data.len().max(1) as f64;
    let bar_w = (slot * 0.58).min(72.0);
    for (index, item) in data.iter().enumerate() {
        let cx = x + slot * (index as f64 + 0.5);
        let bar_h = item.ratio.clamp(0.0, 100.0) / 100.0 * height;
        let color = if index == 0 { SLATE } else { BLUE };
        rect(out, cx - bar_w / 2.0, baseline - bar_h, bar_w, bar_h, color, "#FFFFFF", 0.8, 0.0);
        text(out, cx, baseline - bar_h - 8.0, 11.0, "middle", INK, "normal", &format!("{:.1} ({:.1}%)", item.total, item.ratio));
        text_rot(out, cx, baseline + 18.0, 11.0, "end", INK, &item.label, 25.0);
    }
    line(out, x, baseline, x + width, baseline, SLATE, 1.0);
    text_rot(out, x - 35.0, y + height / 2.0, 12.0, "middle", INK, tr(lang, "总分一致率 (%)", "Score ratio (%)"), -90.0);
    for value in [0, 50, 100] {
        let py = baseline - value as f64 / 100.0 * height;
        line(out, x - 6.0, py, x, py, SLATE, 1.0);
        text(out, x - 10.0, py + 4.0, 11.0, "end", INK, "normal", &value.to_string());
    }
    let reference_label = if lang == "zh" { format!("参考：{}", reference) } else { format!("Reference: {}", reference) };
    text(out, x + width, y + 10.0, 11.0, "end", SLATE, "normal", &reference_label);
}

fn mtba_svg(paralogs: &[Composition], orthologs: &[Composition], lang: &str) -> String {
    let width = 1500.0;
    let height = 1060.0;
    let mut svg = String::new();
    let _ = writeln!(
        svg,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{:.0}\" height=\"{:.0}\" viewBox=\"0 0 {:.0} {:.0}\">",
        width, height, width, height
    );
    rect(&mut svg, 0.0, 0.0, width, height, "#FFFFFF", "#FFFFFF", 0.0, 0.0);
    composition_panel(&mut svg, 110.0, 120.0, 560.0, 300.0, tr(lang, "A  旁系同源：二级结构组成", "A  Paralogs: secondary-structure composition"), paralogs, lang);
    composition_panel(&mut svg, 820.0, 120.0, 560.0, 300.0, tr(lang, "B  直系同源：二级结构组成", "B  Orthologs: secondary-structure composition"), orthologs, lang);
    score_panel(&mut svg, 110.0, 610.0, 560.0, 300.0, tr(lang, "C  旁系同源：总分一致率", "C  Paralogs: score ratio"), paralogs, "MtqA", lang);
    score_panel(&mut svg, 820.0, 610.0, 560.0, 300.0, tr(lang, "D  直系同源：总分一致率", "D  Orthologs: score ratio"), orthologs, "WP_048123177.1", lang);
    let legend_y = 1000.0;
    let legend = [
        (PINK, tr(lang, "α-螺旋 (H)", "α-helix (H)")),
        (BLUE, tr(lang, "β-折叠 (E)", "β-sheet (E)")),
        (GRAY, tr(lang, "coil/卷曲 (C)", "Coil (C)")),
    ];
    let mut x = 520.0;
    for (color, label) in legend {
        rect(&mut svg, x, legend_y - 7.0, 14.0, 14.0, color, "#FFFFFF", 0.5, 0.0);
        text(&mut svg, x + 20.0, legend_y + 4.0, 14.0, "start", INK, "normal", label);
        x += 190.0;
    }
    svg.push_str("</svg>\n");
    svg
}

fn workflow_box(out: &mut String, x: f64, y: f64, w: f64, h: f64, lines: &[&str], fill: &str, font: f64) {
    rect(out, x, y, w, h, fill, SLATE, 2.0, 12.0);
    let base = y + h / 2.0 - (lines.len() as f64 - 1.0) * font * 0.6;
    for (index, line) in lines.iter().enumerate() {
        text(out, x + w / 2.0, base + index as f64 * font * 1.25, font, "middle", INK, "normal", line);
    }
}

fn arrow(out: &mut String, x1: f64, y1: f64, x2: f64, y2: f64) {
    let _ = write!(
        out,
        "<line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"{}\" stroke-width=\"2\" marker-end=\"url(#arrow)\"/>\n",
        x1, y1, x2, y2, SLATE
    );
}

fn workflow_node(
    out: &mut String,
    x: f64,
    y: f64,
    w: f64,
    h: f64,
    title: &str,
    lines: &[&str],
    fill: &str,
) {
    let _ = write!(
        out,
        "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" fill=\"{}\" stroke=\"{}\" stroke-width=\"1.2\"/>\n",
        x, y, w, h, fill, SLATE
    );
    text(out, x + w / 2.0, y + 52.0, 18.0, "middle", INK, "bold", title);
    for (index, value) in lines.iter().enumerate() {
        text(out, x + w / 2.0, y + 92.0 + index as f64 * 28.0, 14.0, "middle", INK, "normal", value);
    }
}

fn workflow_svg(lang: &str) -> String {
    let mut svg = String::new();
    let _ = writeln!(
        svg,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1800\" height=\"800\" viewBox=\"0 0 1800 800\">"
    );
    rect(&mut svg, 0.0, 0.0, 1800.0, 800.0, "#FFFFFF", "#FFFFFF", 0.0, 0.0);
    svg.push_str(r##"<defs><marker id="arrow" markerWidth="9" markerHeight="9" refX="8" refY="4.5" orient="auto"><path d="M0,0 L9,4.5 L0,9 z" fill="#94A3B8"/></marker></defs>"##);
    text(&mut svg, 900.0, 58.0, 22.0, "middle", INK, "bold", tr(lang, "protstruct–PSIPRED 工作流", "protstruct–PSIPRED workflow"));

    workflow_node(&mut svg, 70.0, 220.0, 250.0, 180.0, "FASTA input", &["1–20 sequences", ".fa / .faa / .fasta"], PALE);
    workflow_node(&mut svg, 390.0, 220.0, 270.0, 180.0, "PSIPRED V4.1", &["Rust inference", "H / E / C"], PINK);
    workflow_node(&mut svg, 710.0, 220.0, 370.0, 180.0, "Feature + alignment", &["Disulfide · salt bridge", "MSA · optional PDB"], MINT);
    workflow_node(&mut svg, 1150.0, 220.0, 250.0, 180.0, "2D scoring", &["SS +1 / mismatch 0", "Feature pair +0.5"], PEACH);
    workflow_node(&mut svg, 1470.0, 220.0, 260.0, 180.0, "Results", &["CSV / report / SVG", "Search TSV"], BLUE);

    arrow(&mut svg, 326.0, 310.0, 382.0, 310.0);
    arrow(&mut svg, 666.0, 310.0, 702.0, 310.0);
    arrow(&mut svg, 1086.0, 310.0, 1142.0, 310.0);
    arrow(&mut svg, 1406.0, 310.0, 1462.0, 310.0);

    let _ = write!(
        svg,
        "<path d=\"M525 402 L525 520 Q525 548 553 548 L820 548\" fill=\"none\" stroke=\"#94A3B8\" stroke-width=\"1.5\" stroke-dasharray=\"5 5\" marker-end=\"url(#arrow)\"/>\n"
    );
    let _ = write!(
        svg,
        "<rect x=\"820\" y=\"500\" width=\"910\" height=\"96\" fill=\"#FFFFFF\" stroke=\"{}\" stroke-width=\"1.2\"/>\n",
        SLATE
    );
    text(&mut svg, 1275.0, 538.0, 16.0, "middle", INK, "bold", "protstruct-search");
    text(&mut svg, 1275.0, 568.0, 13.0, "middle", SLATE, "normal", "SS k-mer prefilter → Smith–Waterman → ranked TSV");
    svg.push_str("</svg>\n");
    svg
}

fn parse_args() -> (String, String, String, String) {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("usage: protstruct-plot <workflow|consistency|mtba> [--input file] [--lang zh|en] --output file.svg");
        process::exit(2);
    }
    let mode = args[1].clone();
    let mut input = String::new();
    let mut output = String::new();
    let mut lang = "en".to_string();
    let mut index = 2;
    while index < args.len() {
        match args[index].as_str() {
            "--input" if index + 1 < args.len() => {
                input = args[index + 1].clone();
                index += 2;
            }
            "--output" if index + 1 < args.len() => {
                output = args[index + 1].clone();
                index += 2;
            }
            "--lang" if index + 1 < args.len() => {
                lang = args[index + 1].clone();
                index += 2;
            }
            _ => index += 1,
        }
    }
    if output.is_empty() {
        eprintln!("missing --output");
        process::exit(2);
    }
    (mode, input, output, lang)
}

fn main() {
    let (mode, input, output, lang) = parse_args();
    let svg = match mode.as_str() {
        "workflow" => workflow_svg(&lang),
        "consistency" => match parse_consistency(&input) {
            Ok(data) => consistency_svg(&data, &lang),
            Err(error) => {
                eprintln!("{}", error);
                process::exit(2);
            }
        },
        "mtba" => match parse_mtba(&input) {
            Ok((paralogs, orthologs)) => mtba_svg(&paralogs, &orthologs, &lang),
            Err(error) => {
                eprintln!("{}", error);
                process::exit(2);
            }
        },
        _ => {
            eprintln!("unknown mode: {}", mode);
            process::exit(2);
        }
    };
    if let Err(error) = fs::write(&output, svg) {
        eprintln!("write {}: {}", output, error);
        process::exit(2);
    }
}
