//! GOR IV/V predictor using the compact gor5_4 + gor5_3 model binaries.
//!
//! Usage:
//!   protstruct-gorv --model-dir ../protstruct/data --input input.fa --output ss.fa
//!
//! The binary model layout is documented in scripts/build_gor_models.py.

use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process;

use rayon::prelude::*;

const AA: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
const WINDOW: usize = 21;
const HALF: usize = 10;
const EPS: f64 = 1e-10;
const M6_SIZE: usize = 3 * 20 * 20 * WINDOW * 20 * WINDOW;
const M4_SIZE: usize = 20 * 3 * 20 * WINDOW;

struct Model {
    m6: Vec<u16>,
    m4: Vec<u16>,
}

#[derive(Clone)]
struct Record {
    id: String,
    sequence: String,
}

fn aa_index(byte: u8) -> Option<usize> {
    AA.iter().position(|residue| *residue == byte)
}

fn load_model(path: &Path, expected: usize) -> Vec<u16> {
    let data = fs::read(path).unwrap_or_else(|error| {
        eprintln!("failed to read {}: {}", path.display(), error);
        process::exit(2);
    });
    if data.len() < 12 || &data[0..4] != b"GORV" {
        eprintln!("bad GOR model header: {}", path.display());
        process::exit(2);
    }
    let version = u16::from_le_bytes([data[4], data[5]]);
    let window = u16::from_le_bytes([data[6], data[7]]);
    let states = u16::from_le_bytes([data[8], data[9]]);
    let amino = u16::from_le_bytes([data[10], data[11]]);
    if version != 1 || window as usize != WINDOW || states != 3 || amino != 20 {
        eprintln!("unexpected GOR model header: {}", path.display());
        process::exit(2);
    }
    let body = &data[12..];
    if body.len() % 2 != 0 || body.len() / 2 != expected {
        eprintln!("unexpected GOR model size: {}", path.display());
        process::exit(2);
    }
    let mut values = Vec::with_capacity(expected);
    for chunk in body.chunks_exact(2) {
        values.push(u16::from_le_bytes([chunk[0], chunk[1]]));
    }
    values
}

fn base6(sec: usize, center: usize, outer: usize, position: usize) -> usize {
    (((sec * 20 + center) * 20 + outer) * WINDOW + position) * (20 * WINDOW)
}

fn base4(center: usize, sec: usize) -> usize {
    ((center * 3 + sec) * 20) * WINDOW
}

fn state_score(window: &[u8], state: usize, model: &Model) -> f64 {
    let center = match aa_index(window[HALF]) {
        Some(value) => value,
        None => return f64::NEG_INFINITY,
    };
    let mut outer_sum = 0.0;
    let mut gor3_sum = 0.0;
    for outer in 0..WINDOW {
        let outer_aa = match aa_index(window[outer]) {
            Some(value) => value,
            None => continue,
        };
        let base = base6(state, center, outer_aa, outer);
        let mut inner_sum = 0.0;
        for inner in outer + 1..WINDOW {
            let inner_aa = match aa_index(window[inner]) {
                Some(value) => value,
                None => continue,
            };
            let offset = inner_aa * WINDOW + inner;
            let count = model.m6[base + offset] as f64;
            let mut other = 0.0;
            for alt in 0..3 {
                if alt != state {
                    other += model.m6[base6(alt, center, outer_aa, outer) + offset] as f64;
                }
            }
            inner_sum += ((count + EPS) / (other + EPS)).ln();
        }
        outer_sum += inner_sum;

        let base4_index = base4(center, state) + outer_aa * WINDOW + outer;
        let count4 = model.m4[base4_index] as f64;
        let mut other4 = 0.0;
        for alt in 0..3 {
            if alt != state {
                other4 += model.m4[base4(center, alt) + outer_aa * WINDOW + outer] as f64;
            }
        }
        gor3_sum += ((count4 + EPS) / (other4 + EPS)).ln();
    }
    outer_sum * (2.0 / WINDOW as f64) - gor3_sum * ((WINDOW as f64 - 2.0) / WINDOW as f64)
}

fn post_process(sequence: &mut [u8]) {
    if sequence.len() < 3 {
        return;
    }
    for index in 1..sequence.len() - 1 {
        if sequence[index] != sequence[index - 1] && sequence[index] != sequence[index + 1] {
            sequence[index] = if sequence[index - 1] == sequence[index + 1] {
                sequence[index - 1]
            } else {
                b'C'
            };
        }
    }
}

fn clean_runs(sequence: &mut [u8], kind: u8, min_length: usize) {
    let mut index = 0;
    while index < sequence.len() {
        if sequence[index] != kind {
            index += 1;
            continue;
        }
        let mut end = index;
        while end < sequence.len() && sequence[end] == kind {
            end += 1;
        }
        if end - index < min_length {
            for position in index..end {
                sequence[position] = b'C';
            }
        }
        index = end;
    }
}

fn predict(sequence: &str, model: &Model, post: bool) -> String {
    let bytes = sequence.as_bytes();
    let mut result = vec![b'C'; bytes.len()];
    if bytes.len() >= WINDOW {
        for center in HALF..bytes.len() - HALF {
            let window = &bytes[center - HALF..center + HALF + 1];
            if aa_index(window[HALF]).is_none() {
                continue;
            }
            let mut best = b'C';
            let mut best_score = f64::NEG_INFINITY;
            for (state_index, state_byte) in [(2usize, b'H'), (1usize, b'E'), (0usize, b'C')].iter() {
                let score = state_score(window, *state_index, model);
                if score > best_score {
                    best_score = score;
                    best = *state_byte;
                }
            }
            result[center] = best;
        }
    }
    if post {
        post_process(&mut result);
    }
    // Proline is not compatible with an alpha-helix interior.
    for (index, byte) in bytes.iter().enumerate() {
        if *byte == b'P' {
            result[index] = b'C';
        }
    }
    clean_runs(&mut result, b'H', 4);
    clean_runs(&mut result, b'E', 3);
    String::from_utf8(result).unwrap_or_else(|_| "C".repeat(bytes.len()))
}

fn read_fasta(path: &Path) -> Vec<Record> {
    let text = fs::read_to_string(path).unwrap_or_else(|error| {
        eprintln!("failed to read {}: {}", path.display(), error);
        process::exit(2);
    });
    let mut records = Vec::new();
    let mut id = String::new();
    let mut sequence = String::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with(';') || line.starts_with('#') {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            if !id.is_empty() {
                records.push(Record { id: id.clone(), sequence: sequence.clone() });
                sequence.clear();
            }
            id = header.split_whitespace().next().unwrap_or("sequence").to_string();
        } else {
            sequence.push_str(&line.chars().filter(|c| c.is_ascii_alphabetic()).collect::<String>());
        }
    }
    if !id.is_empty() {
        records.push(Record { id, sequence });
    }
    records
}

fn parse_args() -> (PathBuf, PathBuf, Option<PathBuf>, bool, usize) {
    let mut model_dir = PathBuf::from("../../protstruct/data");
    let mut input = PathBuf::new();
    let mut output = None;
    let mut post = true;
    let mut threads = 0usize;
    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--model-dir" => model_dir = PathBuf::from(args.next().unwrap_or_default()),
            "--input" => input = PathBuf::from(args.next().unwrap_or_default()),
            "--output" => output = Some(PathBuf::from(args.next().unwrap_or_default())),
            "--no-post" => post = false,
            "--threads" => threads = args.next().and_then(|value| value.parse().ok()).unwrap_or(0),
            _ => {}
        }
    }
    if input.as_os_str().is_empty() {
        eprintln!("usage: protstruct-gorv --input input.fa [--model-dir DIR] [--output out.fa] [--no-post] [--threads N]");
        process::exit(2);
    }
    (model_dir, input, output, post, threads)
}

fn main() {
    let (model_dir, input, output, post, threads) = parse_args();
    let m6 = load_model(&model_dir.join("gor5_4_21.bin"), M6_SIZE);
    let m4 = load_model(&model_dir.join("gor5_3_21.bin"), M4_SIZE);
    let model = Model { m6, m4 };
    let records = read_fasta(&input);

    if threads > 0 {
        if let Err(error) = rayon::ThreadPoolBuilder::new().num_threads(threads).build_global() {
            eprintln!("warning: failed to set threads: {}", error);
        }
    }

    let results: Vec<String> = records
        .par_iter()
        .map(|record| {
            let ss = predict(&record.sequence, &model, post);
            format!(">{}\n{}\n", record.id, ss)
        })
        .collect();
    let text = results.concat();

    if let Some(path) = output {
        if let Err(error) = fs::write(&path, text) {
            eprintln!("failed to write {}: {}", path.display(), error);
            process::exit(2);
        }
    } else {
        print!("{}", text);
    }
}
