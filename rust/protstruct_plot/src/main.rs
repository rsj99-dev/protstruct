//! Dependency-free SVG figures for protstruct.
//!
//! Modes:
//!   protstruct-plot workflow --output fig1.svg
//!   protstruct-plot consistency --input consistency.tsv --output fig2.svg
//!   protstruct-plot mtba --input mtba.tsv --output fig3.svg
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

fn fig2(data: &Consistency) -> String {
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
    text(&mut svg, 800.0, 30.0, 24.0, "middle", INK, "bold", &format!("多序列二级结构一致性图（参考序列：{}）", data.reference));
    text(&mut svg, 800.0, 56.0, 13.0, "middle", SLATE, "normal", "逐列比较：二级结构相同 +1，不同 0；二硫键/盐桥对应 +0.5；参考序列为第一条序列");

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
    text_rot(&mut svg, left - 48.0, agreement_y + agreement_h * 0.55, 12.0, "middle", INK, "逐列一致性", -90.0);

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
            format!("{}（参考）", data.rows[row_index].id)
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
    text(&mut svg, left + map_w / 2.0, map_y + map_h + 42.0, 12.0, "middle", INK, "normal", "比对位置（列）");

    // Score panel.
    let score_max = data.max_total.max(1.0);
    line(&mut svg, score_x, map_y + map_h, score_x + score_w, map_y + map_h, SLATE, 1.0);
    line(&mut svg, score_x, map_y, score_x, map_y + map_h, SLATE, 1.0);
    for row_index in 0..rows {
        let bar_w = (data.rows[row_index].total / score_max).clamp(0.0, 1.0) * score_w;
        let color = if data.rows[row_index].reference { SLATE } else { BLUE };
        rect(&mut svg, score_x, map_y + row_index as f64 * row_h + 12.0, bar_w, row_h - 24.0, color, color, 0.0, 0.0);
        text(&mut svg, score_x + bar_w + 7.0, map_y + (row_index as f64 + 0.5) * row_h + 4.0, 12.0, "start", INK, "normal", &format!("{:.1}", data.rows[row_index].total));
    }
    text(&mut svg, score_x + score_w / 2.0, map_y + map_h + 42.0, 12.0, "middle", INK, "normal", &format!("总分（满分 {:.1}）", data.max_total));

    // Legend.
    let legend_y = map_y + map_h + 78.0;
    let legend = [
        (PINK, "α-螺旋 (H)", "rect"),
        (BLUE, "β-折叠 (E)", "rect"),
        (GRAY, "coil/卷曲 (C)", "rect"),
        (SLATE, "二硫键残基", "circle"),
        (PEACH, "盐桥残基", "square"),
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
        if fields.len() < 6 {
            continue;
        }
        let item = Composition {
            label: fields[1].to_string(),
            h: fields[2].parse().unwrap_or(0.0),
            e: fields[3].parse().unwrap_or(0.0),
            c: fields[4].parse().unwrap_or(0.0),
            total: fields[5].parse().unwrap_or(0.0),
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
    text_rot(out, x - 35.0, y + height / 2.0, 12.0, "middle", INK, "占比 (%)", -90.0);
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
) {
    text(out, x, y - 12.0, 18.0, "start", INK, "bold", title);
    let baseline = y + height;
    let slot = width / data.len().max(1) as f64;
    let bar_w = (slot * 0.58).min(72.0);
    let max_score = data.iter().map(|item| item.total).fold(1.0, f64::max);
    for (index, item) in data.iter().enumerate() {
        let cx = x + slot * (index as f64 + 0.5);
        let bar_h = item.total / max_score * height;
        let color = if index == 0 { SLATE } else { BLUE };
        rect(out, cx - bar_w / 2.0, baseline - bar_h, bar_w, bar_h, color, "#FFFFFF", 0.8, 0.0);
        text(out, cx, baseline - bar_h - 8.0, 12.0, "middle", INK, "normal", &format!("{:.1}", item.total));
        text_rot(out, cx, baseline + 18.0, 11.0, "end", INK, &item.label, 25.0);
    }
    line(out, x, baseline, x + width, baseline, SLATE, 1.0);
    text_rot(out, x - 35.0, y + height / 2.0, 12.0, "middle", INK, "总分", -90.0);
    for value in [0, (max_score / 2.0) as i32, max_score as i32] {
        let py = baseline - value as f64 / max_score * height;
        line(out, x - 6.0, py, x, py, SLATE, 1.0);
        text(out, x - 10.0, py + 4.0, 11.0, "end", INK, "normal", &value.to_string());
    }
    text(out, x + width, y + 10.0, 11.0, "end", SLATE, "normal", &format!("参考：{}", reference));
}

fn fig3(paralogs: &[Composition], orthologs: &[Composition]) -> String {
    let width = 1500.0;
    let height = 1060.0;
    let mut svg = String::new();
    let _ = writeln!(
        svg,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{:.0}\" height=\"{:.0}\" viewBox=\"0 0 {:.0} {:.0}\">",
        width, height, width, height
    );
    rect(&mut svg, 0.0, 0.0, width, height, "#FFFFFF", "#FFFFFF", 0.0, 0.0);
    composition_panel(&mut svg, 110.0, 120.0, 560.0, 300.0, "A  旁系同源：二级结构组成", paralogs);
    composition_panel(&mut svg, 820.0, 120.0, 560.0, 300.0, "B  直系同源：二级结构组成", orthologs);
    score_panel(&mut svg, 110.0, 610.0, 560.0, 300.0, "C  旁系同源：一致性总分", paralogs, "MtqA");
    score_panel(&mut svg, 820.0, 610.0, 560.0, 300.0, "D  直系同源：一致性总分", orthologs, "WP_048123177.1");
    let legend_y = 1000.0;
    let legend = [(PINK, "α-螺旋 (H)"), (BLUE, "β-折叠 (E)"), (GRAY, "coil/卷曲 (C)")];
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

fn fig1() -> String {
    let mut svg = String::new();
    let _ = writeln!(
        svg,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1200\" height=\"760\" viewBox=\"0 0 1200 760\">"
    );
    rect(&mut svg, 0.0, 0.0, 1200.0, 760.0, "#FFFFFF", "#FFFFFF", 0.0, 0.0);
    svg.push_str("<defs><marker id=\"arrow\" markerWidth=\"8\" markerHeight=\"8\" refX=\"7\" refY=\"4\" orient=\"auto\"><path d=\"M0,0 L8,4 L0,8 z\" fill=\"#5878A6\"/></marker></defs>\n");
    workflow_box(&mut svg, 300.0, 40.0, 600.0, 65.0, &["输入 FASTA 文件", ".fa / .faa / .fasta，1–20 条序列"], PALE, 17.0);
    workflow_box(&mut svg, 300.0, 145.0, 600.0, 65.0, &["解析与校验", "编码识别 · 残基校验 · 重复 ID 处理"], GRAY, 17.0);
    rect(&mut svg, 70.0, 270.0, 1060.0, 250.0, BG, GRAY, 1.0, 15.0);
    text(&mut svg, 90.0, 300.0, 14.0, "start", SLATE, "normal", "核心分析模块（Rust/Python，无需联网）");
    workflow_box(&mut svg, 100.0, 330.0, 450.0, 75.0, &["二级结构预测", "GOR IV → H / E / C"], PINK, 17.0);
    workflow_box(&mut svg, 650.0, 330.0, 450.0, 75.0, &["多序列比对", "center-star + NW（BLOSUM62）"], BLUE, 17.0);
    workflow_box(&mut svg, 100.0, 425.0, 450.0, 70.0, &["二硫键预测", "动态规划最优配对 · 游离 Cys 输出"], PEACH, 17.0);
    workflow_box(&mut svg, 650.0, 425.0, 450.0, 70.0, &["盐桥候选", "酸性–碱性残基 · 二级结构上下文打分"], MINT, 17.0);
    workflow_box(&mut svg, 300.0, 565.0, 600.0, 65.0, &["一致性评分", "二级结构相同 +1 · 不同 0 · 二硫键/盐桥对应 +0.5"], GRAY, 16.0);
    workflow_box(&mut svg, 300.0, 665.0, 600.0, 65.0, &["结果输出", "报告 / CSV / SVG / 检索"], PALE, 17.0);
    arrow(&mut svg, 600.0, 105.0, 600.0, 142.0);
    arrow(&mut svg, 600.0, 210.0, 600.0, 267.0);
    arrow(&mut svg, 600.0, 520.0, 600.0, 562.0);
    arrow(&mut svg, 600.0, 630.0, 600.0, 662.0);
    svg.push_str("</svg>\n");
    svg
}

fn parse_args() -> (String, String, String) {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("usage: protstruct-plot <workflow|consistency|mtba> [--input file] --output file.svg");
        process::exit(2);
    }
    let mode = args[1].clone();
    let mut input = String::new();
    let mut output = String::new();
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
            _ => index += 1,
        }
    }
    if output.is_empty() {
        eprintln!("missing --output");
        process::exit(2);
    }
    (mode, input, output)
}

fn main() {
    let (mode, input, output) = parse_args();
    let svg = match mode.as_str() {
        "workflow" => fig1(),
        "consistency" => match parse_consistency(&input) {
            Ok(data) => fig2(&data),
            Err(error) => {
                eprintln!("{}", error);
                process::exit(2);
            }
        },
        "mtba" => match parse_mtba(&input) {
            Ok((paralogs, orthologs)) => fig3(&paralogs, &orthologs),
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
