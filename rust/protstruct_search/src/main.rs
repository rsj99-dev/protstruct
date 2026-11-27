//! 二级结构数据库检索：批量 PSIPRED、k-mer 预筛和并行 Smith–Waterman。

use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::env;
use std::fs;
use std::process;

use protstruct_psipred::PsipredModel;
use rayon::prelude::*;

#[derive(Clone)]
struct Record {
    id: String,
    description: String,
    sequence: String,
    ss: String,
}

#[derive(Clone)]
struct Hit {
    id: String,
    description: String,
    length: usize,
    ss_score: i32,
    matches: i32,
    aligned: i32,
    coverage: f64,
    ss_identity: f64,
    seq_sim: f64,
    total_score: f64,
}

fn fail(message: &str) -> ! {
    eprintln!("{}", message);
    process::exit(2);
}

fn clean_sequence(raw: &str) -> String {
    let mut cleaned = String::new();
    for ch in raw.chars() {
        if ch.is_whitespace() || ch == '-' || ch == '.' || ch == '*' {
            continue;
        }
        let upper = ch.to_ascii_uppercase();
        let mapped = match upper {
            'U' => 'C',
            'O' => 'K',
            'B' | 'Z' | 'J' | 'X' => 'X',
            value if "ACDEFGHIKLMNPQRSTVWY".contains(value) => value,
            _ => fail(&format!("invalid residue: {}", ch)),
        };
        cleaned.push(mapped);
    }
    cleaned
}

fn parse_fasta(path: &str) -> Vec<Record> {
    let text = fs::read_to_string(path).unwrap_or_else(|error| fail(&format!("read {}: {}", path, error)));
    let mut records: Vec<Record> = Vec::new();
    let mut seen: HashMap<String, usize> = HashMap::new();
    let mut header: Option<String> = None;
    let mut sequence = String::new();
    for raw_line in text.trim_start_matches('\u{feff}').lines() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with(';') || line.starts_with('#') {
            continue;
        }
        if let Some(value) = line.strip_prefix('>') {
            if let Some(previous) = header.take() {
                let (id, description) = split_header(&previous);
                records.push(Record { id, description, sequence: clean_sequence(&sequence), ss: String::new() });
                sequence.clear();
            }
            let value = value.trim();
            if value.is_empty() {
                fail("empty FASTA header");
            }
            header = Some(value.to_string());
        } else if header.is_some() {
            sequence.push_str(line);
        } else {
            fail("sequence data found before the first '>' header");
        }
    }
    if let Some(previous) = header {
        let (id, description) = split_header(&previous);
        records.push(Record { id, description, sequence: clean_sequence(&sequence), ss: String::new() });
    }
    if records.is_empty() {
        fail("no sequences found");
    }
    for record in &mut records {
        if record.sequence.is_empty() {
            fail(&format!("sequence {} is empty", record.id));
        }
        let count = seen.entry(record.id.clone()).or_insert(0);
        *count += 1;
        if *count > 1 {
            record.id = format!("{}#{}", record.id, count);
        }
    }
    records
}

fn split_header(header: &str) -> (String, String) {
    let mut fields = header.splitn(2, char::is_whitespace);
    let id = fields.next().unwrap_or("sequence").to_string();
    let description = fields.next().unwrap_or("").trim().to_string();
    (id, description)
}

fn kmer_set(sequence: &str, k: usize) -> HashSet<String> {
    if k == 0 || sequence.len() < k {
        return HashSet::new();
    }
    (0..=sequence.len() - k).map(|index| sequence[index..index + k].to_string()).collect()
}

fn jaccard(first: &HashSet<String>, second: &HashSet<String>) -> f64 {
    if first.is_empty() || second.is_empty() {
        return 0.0;
    }
    let intersection = first.intersection(second).count() as f64;
    let union = first.union(second).count() as f64;
    if union == 0.0 { 0.0 } else { intersection / union }
}

/// H/E/C 字母表上的 Smith–Waterman：match +2，mismatch -1，gap -2。
fn smith_waterman(query: &str, target: &str) -> (i32, i32, i32) {
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
            let substitution = if query.as_bytes()[i - 1] == target.as_bytes()[j - 1] { 2 } else { -1 };
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
        let substitution = if query.as_bytes()[i - 1] == target.as_bytes()[j - 1] { 2 } else { -1 };
        if current == matrix[i - 1][j - 1] + substitution {
            aligned += 1;
            if query.as_bytes()[i - 1] == target.as_bytes()[j - 1] {
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

#[derive(Default)]
struct Options {
    input: String,
    query: String,
    database: String,
    data: String,
    top: usize,
    prefilter: usize,
    kmer: usize,
    min_ss_identity: f64,
    min_coverage: f64,
}

fn parse_args() -> Options {
    let mut options = Options { top: 20, prefilter: 2000, kmer: 3, ..Options::default() };
    let mut args = env::args().skip(1);
    while let Some(arg) = args.next() {
        let mut value = || args.next().unwrap_or_default();
        match arg.as_str() {
            "--input" => options.input = value(),
            "--query" => options.query = value(),
            "--database" => options.database = value(),
            "--data" => options.data = value(),
            "--top" => options.top = value().parse().unwrap_or(20),
            "--prefilter" => options.prefilter = value().parse().unwrap_or(2000),
            "--kmer" => options.kmer = value().parse().unwrap_or(3),
            "--min-ss-identity" => options.min_ss_identity = value().parse().unwrap_or(0.0),
            "--min-coverage" => options.min_coverage = value().parse().unwrap_or(0.0),
            "--help" | "-h" => {
                println!("usage: protstruct-search --query <query.fa> --database <database.fa> --data <weights-dir> [--top N] [--prefilter N] [--kmer K]");
                process::exit(0);
            }
            _ => fail(&format!("unknown argument: {}", arg)),
        }
    }
    if options.top == 0 || options.prefilter == 0 || options.kmer == 0 {
        fail("top, prefilter, and kmer must be positive");
    }
    options
}

fn legacy_score(path: &str, top: usize) {
    let content = fs::read_to_string(path).unwrap_or_else(|error| fail(&format!("read {}: {}", path, error)));
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
            description: String::new(),
            ss: fields[1].to_string(),
            sequence: fields[2].to_string(),
        };
        if record.id == "__QUERY__" {
            query = Some(record);
        } else {
            records.push(record);
        }
    }
    let query = query.unwrap_or_else(|| fail("missing query record"));
    let query_kmers = kmer_set(&query.sequence, 3);
    let mut hits: Vec<(String, i32, i32, i32, f64)> = records
        .par_iter()
        .map(|record| {
            let (score, matches, aligned) = smith_waterman(&query.ss, &record.ss);
            let seq_sim = jaccard(&query_kmers, &kmer_set(&record.sequence, 3));
            (record.id.clone(), score, matches, aligned, seq_sim)
        })
        .collect();
    hits.sort_by(|left, right| {
        let left_score = left.1 as f64 + 2.0 * left.4;
        let right_score = right.1 as f64 + 2.0 * right.4;
        right_score.partial_cmp(&left_score).unwrap_or(Ordering::Equal)
    });
    hits.truncate(top);
    for hit in hits {
        println!("{}\t{}\t{}\t{}\t{:.4}", hit.0, hit.1, hit.2, hit.3, hit.4);
    }
}

fn postprocess_states(sequence: &str, states: &str) -> String {
    let mut chars: Vec<u8> = states.as_bytes().to_vec();
    for (index, residue) in sequence.bytes().enumerate() {
        if residue == b'P' {
            chars[index] = b'C';
        }
    }
    clean_runs(&mut chars, b'H', 4);
    clean_runs(&mut chars, b'E', 3);
    String::from_utf8(chars).unwrap_or_else(|_| fail("internal state encoding error"))
}

fn clean_runs(chars: &mut [u8], kind: u8, min_length: usize) {
    let mut index = 0;
    while index < chars.len() {
        if chars[index] != kind {
            index += 1;
            continue;
        }
        let end = (index + 1..chars.len()).find(|position| chars[*position] != kind).unwrap_or(chars.len());
        if end - index < min_length {
            for value in &mut chars[index..end] {
                *value = b'C';
            }
        }
        index = end;
    }
}

fn full_search(options: &Options) {
    if options.query.is_empty() || options.database.is_empty() || options.data.is_empty() {
        fail("full search requires --query, --database, and --data");
    }
    let model = PsipredModel::load(&options.data).unwrap_or_else(|error| fail(&error));
    let mut database = parse_fasta(&options.database);
    let mut queries = parse_fasta(&options.query);
    let query = queries.remove(0);
    let query_ss = postprocess_states(
        &query.sequence,
        &model.predict(&query.sequence).unwrap_or_else(|error| fail(&error)),
    );

    database.par_iter_mut().for_each(|record| {
        let states = model.predict(&record.sequence).unwrap_or_else(|error| fail(&error));
        record.ss = postprocess_states(&record.sequence, &states);
    });

    let query_seq_set = kmer_set(&query.sequence, options.kmer);
    let query_ss_set = kmer_set(&query_ss, options.kmer);
    let mut entries: Vec<(f64, usize)> = database
        .iter()
        .enumerate()
        .map(|(index, record)| {
            let seq_similarity = jaccard(&query_seq_set, &kmer_set(&record.sequence, options.kmer));
            let ss_similarity = jaccard(&query_ss_set, &kmer_set(&record.ss, options.kmer));
            (0.55 * ss_similarity + 0.45 * seq_similarity, index)
        })
        .collect();
    entries.sort_by(|left, right| right.0.partial_cmp(&left.0).unwrap_or(Ordering::Equal));
    let candidate_count = options.prefilter.min(entries.len());
    let candidates: Vec<&Record> = entries[..candidate_count].iter().map(|(_, index)| &database[*index]).collect();
    let query_ss_len = query_ss.len();

    let mut hits: Vec<Hit> = candidates
        .par_iter()
        .filter_map(|record| {
            let (score, matches, aligned) = smith_waterman(&query_ss, &record.ss);
            let coverage = aligned as f64 / query_ss_len as f64;
            let identity = if aligned == 0 { 0.0 } else { matches as f64 / aligned as f64 };
            if identity < options.min_ss_identity || coverage < options.min_coverage {
                return None;
            }
            let seq_sim = jaccard(&query_seq_set, &kmer_set(&record.sequence, options.kmer));
            Some(Hit {
                id: record.id.clone(),
                description: record.description.clone(),
                length: record.sequence.len(),
                ss_score: score,
                matches,
                aligned,
                coverage,
                ss_identity: identity,
                seq_sim,
                total_score: score as f64 + 2.0 * seq_sim,
            })
        })
        .collect();
    hits.sort_by(|left, right| {
        right.total_score
            .partial_cmp(&left.total_score)
            .unwrap_or(Ordering::Equal)
            .then(right.ss_identity.partial_cmp(&left.ss_identity).unwrap_or(Ordering::Equal))
    });
    hits.truncate(options.top);

    println!("rank\tidentifier\ttotal_score\tss_score\tss_identity\tcoverage\tseq_similarity\tlength\tdescription\tmatches\taligned_length");
    for (index, hit) in hits.iter().enumerate() {
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            index + 1, hit.id, hit.total_score, hit.ss_score, hit.ss_identity,
            hit.coverage, hit.seq_sim, hit.length, hit.description.replace('\t', " "),
            hit.matches, hit.aligned
        );
    }
}

fn main() {
    let options = parse_args();
    if options.input.is_empty() {
        full_search(&options);
    } else {
        legacy_score(&options.input, options.top);
    }
}
