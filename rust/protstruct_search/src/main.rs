//! protstruct-search: parallel secondary-structure search core.
//!
//! Input (tab-separated):
//!   __QUERY__ <query_ss> <query_sequence>
//!   <id> <target_ss> <target_sequence>
//!
//! Output (tab-separated, ranked):
//!   id  ss_score  matches  aligned_length  sequence_kmer_similarity

use std::collections::HashSet;
use std::env;
use std::fs;
use std::process;

use rayon::prelude::*;

#[derive(Clone)]
struct Record {
    id: String,
    ss: String,
    sequence: String,
}

#[derive(Clone)]
struct Hit {
    id: String,
    ss_score: i32,
    matches: i32,
    aligned: i32,
    seq_sim: f64,
}

fn kmer_set(sequence: &str, k: usize) -> HashSet<String> {
    let chars: Vec<char> = sequence.chars().collect();
    let mut set = HashSet::new();
    if chars.len() < k {
        return set;
    }
    for index in 0..=chars.len() - k {
        set.insert(chars[index..index + k].iter().collect());
    }
    set
}

fn jaccard(first: &HashSet<String>, second: &HashSet<String>) -> f64 {
    if first.is_empty() || second.is_empty() {
        return 0.0;
    }
    let intersection = first.intersection(second).count() as f64;
    let union = first.union(second).count() as f64;
    if union == 0.0 { 0.0 } else { intersection / union }
}

/// Smith-Waterman on the H/E/C alphabet: match +2, mismatch -1, gap -2.
fn smith_waterman(query: &[u8], target: &[u8]) -> (i32, i32, i32) {
    let rows = query.len();
    let cols = target.len();
    if rows == 0 || cols == 0 {
        return (0, 0, 0);
    }
    let mut matrix = vec![vec![0i32; cols + 1]; rows + 1];
    let mut best = 0i32;
    let mut best_i = 0usize;
    let mut best_j = 0usize;

    for i in 1..=rows {
        for j in 1..=cols {
            let substitution = if query[i - 1] == target[j - 1] { 2 } else { -1 };
            let value = 0
                .max(matrix[i - 1][j - 1] + substitution)
                .max(matrix[i - 1][j] - 2)
                .max(matrix[i][j - 1] - 2);
            matrix[i][j] = value;
            if value > best {
                best = value;
                best_i = i;
                best_j = j;
            }
        }
    }

    let mut i = best_i;
    let mut j = best_j;
    let mut matches = 0i32;
    let mut aligned = 0i32;
    while i > 0 && j > 0 && matrix[i][j] > 0 {
        let current = matrix[i][j];
        let substitution = if query[i - 1] == target[j - 1] { 2 } else { -1 };
        if current == matrix[i - 1][j - 1] + substitution {
            aligned += 1;
            if query[i - 1] == target[j - 1] {
                matches += 1;
            }
            i -= 1;
            j -= 1;
        } else if current == matrix[i - 1][j] - 2 {
            aligned += 1;
            i -= 1;
        } else if current == matrix[i][j - 1] - 2 {
            aligned += 1;
            j -= 1;
        } else {
            break;
        }
    }
    (best, matches, aligned)
}

fn score_record(query: &Record, query_kmers: &HashSet<String>, record: &Record) -> Hit {
    let (score, matches, aligned) = smith_waterman(query.ss.as_bytes(), record.ss.as_bytes());
    let seq_sim = jaccard(query_kmers, &kmer_set(&record.sequence, 3));
    Hit {
        id: record.id.clone(),
        ss_score: score,
        matches,
        aligned,
        seq_sim,
    }
}

fn parse_args() -> (String, usize) {
    let mut input = String::new();
    let mut top = 20usize;
    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--input" => input = args.next().unwrap_or_default(),
            "--top" => top = args.next().and_then(|value| value.parse().ok()).unwrap_or(20),
            _ => {}
        }
    }
    if input.is_empty() {
        eprintln!("usage: protstruct-search --input <tsv> [--top N]");
        process::exit(2);
    }
    (input, top)
}

fn main() {
    let (path, top) = parse_args();
    let content = match fs::read_to_string(&path) {
        Ok(value) => value,
        Err(error) => {
            eprintln!("failed to read {}: {}", path, error);
            process::exit(2);
        }
    };

    let mut query: Option<Record> = None;
    let mut records: Vec<Record> = Vec::new();
    for line in content.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            continue;
        }
        let record = Record {
            id: fields[0].to_string(),
            ss: fields[1].to_string(),
            sequence: fields[2].to_string(),
        };
        if record.id == "__QUERY__" {
            query = Some(record);
        } else {
            records.push(record);
        }
    }

    let query = match query {
        Some(value) => value,
        None => {
            eprintln!("missing query record");
            process::exit(2);
        }
    };

    let query_kmers = kmer_set(&query.sequence, 3);
    let mut hits: Vec<Hit> = records
        .par_iter()
        .map(|record| score_record(&query, &query_kmers, record))
        .collect();

    hits.sort_by(|left, right| {
        let left_score = left.ss_score as f64 + 2.0 * left.seq_sim;
        let right_score = right.ss_score as f64 + 2.0 * right.seq_sim;
        right_score
            .partial_cmp(&left_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    hits.truncate(top);

    for hit in hits {
        println!(
            "{}\t{}\t{}\t{}\t{:.4}",
            hit.id, hit.ss_score, hit.matches, hit.aligned, hit.seq_sim
        );
    }
}
