//! 批量两两比对与二硫键配对核心。

use std::cmp::Ordering;
use std::io::{self, Read};
use std::process;
use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};
use std::sync::Mutex;

const NEGATIVE: f64 = -1e9;
const MAX_THREADS: usize = 64;
const STATE_M: u8 = 0;
const STATE_X: u8 = 1;
const STATE_Y: u8 = 2;

// 与 Python API 完全一致的标准 BLOSUM62（C-C = 9）。
const BLOSUM62: [[i8; 20]; 20] = [
    [4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0],
    [-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3],
    [-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3],
    [-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3],
    [0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],
    [-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2],
    [-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2],
    [0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3],
    [-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3],
    [-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3],
    [-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1],
    [-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2],
    [-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1],
    [-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1],
    [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2],
    [1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2],
    [0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0],
    [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3],
    [-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1],
    [0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4],
];

fn fail(message: &str) -> ! {
    eprintln!("{}", message);
    process::exit(2);
}

/// 固定工作线程数、保持结果顺序的轻量任务池。
fn parallel_map<T, R, F>(items: &[T], threads: usize, job: F) -> Vec<R>
where
    T: Sync,
    R: Send,
    F: Fn(&T) -> R + Sync,
{
    if items.is_empty() {
        return Vec::new();
    }
    let workers = threads.max(1).min(MAX_THREADS).min(items.len());
    if workers <= 1 {
        return items.iter().map(job).collect();
    }
    let next = AtomicUsize::new(0);
    let results = Mutex::new(Vec::with_capacity(items.len()));
    std::thread::scope(|scope| {
        for _ in 0..workers {
            scope.spawn(|| loop {
                let index = next.fetch_add(1, AtomicOrdering::Relaxed);
                if index >= items.len() {
                    break;
                }
                let value = job(&items[index]);
                results.lock().unwrap().push((index, value));
            });
        }
    });
    let mut collected = results.into_inner().unwrap();
    collected.sort_by_key(|(index, _)| *index);
    collected.into_iter().map(|(_, value)| value).collect()
}

fn aa_index(ch: u8) -> Option<usize> {
    let upper = ch.to_ascii_uppercase();
    const ORDER: &[u8] = b"ARNDCQEGHILKMFPSTWYV";
    ORDER.iter().position(|value| *value == upper)
}

fn substitution(left: u8, right: u8) -> f64 {
    match (aa_index(left), aa_index(right)) {
        (Some(first), Some(second)) => BLOSUM62[first][second] as f64,
        _ => -1.0,
    }
}

fn best_state(first: f64, second: f64, third: f64) -> u8 {
    let mut best = first;
    let mut state = STATE_M;
    if second > best {
        best = second;
        state = STATE_X;
    }
    if third > best {
        state = STATE_Y;
    }
    state
}

/// 与 Python 版保持相同回退顺序的仿射空位全局比对。
fn align_pair(first: &str, second: &str, gap_open: f64, gap_extend: f64) -> (String, String, f64) {
    let rows = first.len();
    let cols = second.len();
    let stride = cols + 1;
    let cells = (rows + 1) * stride;
    let mut previous_m = vec![NEGATIVE; stride];
    let mut previous_x = vec![NEGATIVE; stride];
    let mut previous_y = vec![NEGATIVE; stride];
    let mut current_m = vec![NEGATIVE; stride];
    let mut current_x = vec![NEGATIVE; stride];
    let mut current_y = vec![NEGATIVE; stride];
    let mut trace_m = vec![STATE_M; cells];
    let mut trace_x = vec![STATE_M; cells];
    let mut trace_y = vec![STATE_M; cells];
    previous_m[0] = 0.0;

    for j in 1..=cols {
        previous_y[j] = gap_open + (j as f64 - 1.0) * gap_extend;
        trace_y[j] = if j == 1 { STATE_M } else { STATE_Y };
    }

    for i in 1..=rows {
        current_m[0] = NEGATIVE;
        current_x[0] = gap_open + (i as f64 - 1.0) * gap_extend;
        current_y[0] = NEGATIVE;
        trace_x[i * stride] = if i == 1 { STATE_M } else { STATE_X };
        for j in 1..=cols {
            let diagonal = best_state(previous_m[j - 1], previous_x[j - 1], previous_y[j - 1]);
            current_m[j] = match diagonal {
                STATE_M => previous_m[j - 1],
                STATE_X => previous_x[j - 1],
                _ => previous_y[j - 1],
            } + substitution(first.as_bytes()[i - 1], second.as_bytes()[j - 1]);
            trace_m[i * stride + j] = diagonal;

            let up = best_state(previous_m[j] + gap_open, previous_x[j] + gap_extend, previous_y[j] + gap_open);
            current_x[j] = match up {
                STATE_M => previous_m[j] + gap_open,
                STATE_X => previous_x[j] + gap_extend,
                _ => previous_y[j] + gap_open,
            };
            trace_x[i * stride + j] = up;

            // Python 版在 Y 状态按 M、Y、X 顺序打破同分。
            let left_values = [
                current_m[j - 1] + gap_open,
                current_y[j - 1] + gap_extend,
                current_x[j - 1] + gap_open,
            ];
            let mut left = STATE_M;
            if left_values[1] > left_values[0] {
                left = STATE_Y;
            }
            let left_best = if left == STATE_M { left_values[0] } else { left_values[1] };
            if left_values[2] > left_best {
                left = STATE_X;
            }
            current_y[j] = match left {
                STATE_M => current_m[j - 1] + gap_open,
                STATE_Y => current_y[j - 1] + gap_extend,
                _ => current_x[j - 1] + gap_open,
            };
            trace_y[i * stride + j] = left;
        }
        std::mem::swap(&mut previous_m, &mut current_m);
        std::mem::swap(&mut previous_x, &mut current_x);
        std::mem::swap(&mut previous_y, &mut current_y);
    }

    let best = previous_m[cols].max(previous_x[cols]).max(previous_y[cols]);
    let mut state = if best == previous_m[cols] {
        STATE_M
    } else if best == previous_x[cols] {
        STATE_X
    } else {
        STATE_Y
    };
    let mut left_chars = Vec::with_capacity(rows + cols);
    let mut right_chars = Vec::with_capacity(rows + cols);
    let mut i = rows;
    let mut j = cols;
    while i > 0 || j > 0 {
        match state {
            STATE_M if i > 0 && j > 0 => {
                left_chars.push(first.as_bytes()[i - 1]);
                right_chars.push(second.as_bytes()[j - 1]);
                state = trace_m[i * stride + j];
                i -= 1;
                j -= 1;
            }
            STATE_X if i > 0 => {
                left_chars.push(first.as_bytes()[i - 1]);
                right_chars.push(b'-');
                state = trace_x[i * stride + j];
                i -= 1;
            }
            _ if j > 0 => {
                left_chars.push(b'-');
                right_chars.push(second.as_bytes()[j - 1]);
                state = trace_y[i * stride + j];
                j -= 1;
            }
            _ => break,
        }
    }
    left_chars.reverse();
    right_chars.reverse();
    (
        String::from_utf8(left_chars).expect("internal alignment encoding error"),
        String::from_utf8(right_chars).expect("internal alignment encoding error"),
        best,
    )
}

fn run_alignment(input: &str, threads: usize) {
    let mut lines = input.lines();
    let header = lines.next().unwrap_or_default();
    let fields: Vec<&str> = header.split('\t').collect();
    if fields.len() != 3 {
        fail("alignment header must be: gap_open gap_extend k");
    }
    let gap_open: f64 = fields[0].parse().unwrap_or_else(|_| fail("invalid gap_open"));
    let gap_extend: f64 = fields[1].parse().unwrap_or_else(|_| fail("invalid gap_extend"));
    let _k: usize = fields[2].parse().unwrap_or_else(|_| fail("invalid k"));
    let records: Vec<(usize, String, String)> = lines
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() != 3 {
                fail("alignment rows must be: index sequence1 sequence2");
            }
            (fields[0].parse().unwrap_or_else(|_| fail("invalid alignment index")), fields[1].to_string(), fields[2].to_string())
        })
        .collect();
    let results: Vec<(usize, String, String, f64)> = parallel_map(&records, threads, |(index, first, second)| {
        let (left, right, score) = align_pair(first, second, gap_open, gap_extend);
        (*index, left, right, score)
    });
    let mut ordered = results;
    ordered.sort_by_key(|record| record.0);
    for (index, left, right, score) in ordered {
        println!("{}\t{}\t{}\t{}", index, left, right, score);
    }
}

fn pair_cost(separation: usize, in_helix: bool, max_loop: usize) -> f64 {
    let mut cost = (separation as f64 - 8.0).abs() / 8.0;
    if separation > max_loop {
        cost += 2.0;
    }
    if in_helix {
        cost += 0.35;
    }
    cost
}

fn solve_pairs(
    mask: u32,
    positions: &[usize],
    costs: &[Option<f64>],
    memo: &mut Vec<Option<(usize, f64)>>,
) -> (usize, f64) {
    if mask == 0 {
        return (0, 0.0);
    }
    if let Some(value) = memo[mask as usize] {
        return value;
    }
    let lowest = mask.trailing_zeros() as usize;
    let rest = mask ^ (1 << lowest);
    let mut best = solve_pairs(rest, positions, costs, memo);
    let mut remaining = rest;
    while remaining != 0 {
        let partner = remaining.trailing_zeros() as usize;
        remaining ^= 1 << partner;
        if let Some(cost) = costs[lowest * positions.len() + partner] {
            let (sub_count, sub_cost) = solve_pairs(rest ^ (1 << partner), positions, costs, memo);
            let candidate = (sub_count + 1, sub_cost + cost);
            if candidate.0 > best.0 || (candidate.0 == best.0 && candidate.1 < best.1 - 1e-12) {
                best = candidate;
            }
        }
    }
    memo[mask as usize] = Some(best);
    best
}

fn optimal_pairs(positions: &[usize], costs: &[Option<f64>], limit: usize) -> Vec<(usize, usize)> {
    let count = positions.len();
    if count > limit {
        let mut candidates: Vec<(f64, usize, usize)> = costs
            .iter()
            .enumerate()
            .filter_map(|(index, cost)| {
                cost.map(|value| (value, index / count, index % count))
            })
            .filter(|(_, first, second)| first < second)
            .collect();
        candidates.sort_by(|left, right| left.0.partial_cmp(&right.0).unwrap_or(Ordering::Equal).then(left.1.cmp(&right.1)).then(left.2.cmp(&right.2)));
        let mut used = vec![false; count];
        let mut pairs = Vec::new();
        for (_, first, second) in candidates {
            if !used[first] && !used[second] {
                used[first] = true;
                used[second] = true;
                pairs.push((first, second));
            }
        }
        return pairs;
    }

    let mut memo = vec![None; 1usize << count];
    let (mut best_count, mut best_cost) = solve_pairs((1 << count) - 1, positions, costs, &mut memo);
    let mut pairs = Vec::new();
    let mut mask = (1u32 << count) - 1;
    while mask != 0 {
        let lowest = mask.trailing_zeros() as usize;
        let rest = mask ^ (1 << lowest);
        let (rest_count, rest_cost) = solve_pairs(rest, positions, costs, &mut memo);
        if rest_count == best_count && (rest_cost - best_cost).abs() < 1e-9 {
            mask = rest;
            continue;
        }
        let mut matched = false;
        let mut remaining = rest;
        while remaining != 0 {
            let partner = remaining.trailing_zeros() as usize;
            remaining ^= 1 << partner;
            if let Some(cost) = costs[lowest * count + partner] {
                let (sub_count, sub_cost) = solve_pairs(rest ^ (1 << partner), positions, costs, &mut memo);
                if sub_count + 1 == best_count && (sub_cost + cost - best_cost).abs() < 1e-9 {
                    pairs.push((lowest, partner));
                    mask = rest ^ (1 << partner);
                    best_count = sub_count;
                    best_cost = sub_cost;
                    matched = true;
                    break;
                }
            }
        }
        if !matched {
            mask = rest;
        }
    }
    pairs
}

fn run_disulfide(input: &str, threads: usize) {
    let records: Vec<(usize, String, String, usize, usize)> = input
        .lines()
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() != 5 {
                fail("disulfide rows must be: index sequence ss min_separation max_loop");
            }
            (
                fields[0].parse().unwrap_or_else(|_| fail("invalid disulfide index")),
                fields[1].to_ascii_uppercase(),
                fields[2].to_string(),
                fields[3].parse().unwrap_or_else(|_| fail("invalid min_separation")),
                fields[4].parse().unwrap_or_else(|_| fail("invalid max_loop")),
            )
        })
        .collect();
    let results: Vec<(usize, Vec<(usize, usize)>)> = parallel_map(&records, threads, |(index, sequence, ss, min_separation, max_loop)| {
            let positions: Vec<usize> = sequence.bytes().enumerate().filter_map(|(position, residue)| if residue == b'C' { Some(position + 1) } else { None }).collect();
            let count = positions.len();
            let mut costs = vec![None; count * count];
            for first in 0..count {
                for second in first + 1..count {
                    let separation = positions[second] - positions[first];
                    if separation < *min_separation {
                        continue;
                    }
                    let in_helix = ss.len() == sequence.len()
                        && ss.as_bytes()[positions[first] - 1] == b'H'
                        && ss.as_bytes()[positions[second] - 1] == b'H';
                    costs[first * count + second] = Some(pair_cost(separation, in_helix, *max_loop));
                }
            }
        let pairs = optimal_pairs(&positions, &costs, 16)
            .into_iter()
            .map(|(first, second)| (positions[first], positions[second]))
            .collect();
        (*index, pairs)
    });
    let mut ordered = results;
    ordered.sort_by_key(|record| record.0);
    for (index, pairs) in ordered {
        if pairs.is_empty() {
            println!("{}", index);
        } else {
            for (first, second) in pairs {
                println!("{}\t{}\t{}", index, first, second);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[ignore]
    fn benchmark_alignment_thread_scaling() {
        let alphabet = b"ACDEFGHIKLMNPQRSTVWY";
        let length = 1500usize;
        let base: Vec<u8> = (0..length).map(|index| alphabet[(index * 7 + index / 13) % 20]).collect();
        let sequences: Vec<String> = (0..16).map(|copy| {
            let mut value = base.clone();
            for position in (copy..length).step_by(19) {
                value[position] = alphabet[(position + copy * 3 + 1) % 20];
            }
            String::from_utf8(value).unwrap()
        }).collect();
        let records: Vec<(usize, String, String)> = (1..sequences.len())
            .map(|index| (index, sequences[0].clone(), sequences[index].clone()))
            .collect();
        for threads in [1usize, 4, 8, 16, 64] {
            let started = std::time::Instant::now();
            let results = parallel_map(&records, threads, |(_, first, second)| align_pair(first, second, -11.0, -1.0));
            println!("threads={:<2} time={:?} pairs={}", threads, started.elapsed(), results.len());
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut mode = String::new();
    let mut requested_threads = 0usize;
    let mut index = 0;
    while index < args.len() {
        match args[index].as_str() {
            "--threads" if index + 1 < args.len() => {
                requested_threads = args[index + 1].parse().unwrap_or_else(|_| fail("invalid --threads value"));
                index += 2;
            }
            "--help" | "-h" => {
                println!("usage: protstruct-core <align|disulfide> [--threads N] < input.tsv");
                return;
            }
            value if mode.is_empty() => {
                mode = value.to_string();
                index += 1;
            }
            _ => fail(&format!("unknown argument: {}", args[index])),
        }
    }
    let available = std::thread::available_parallelism().map(|value| value.get()).unwrap_or(1);
    let threads = if requested_threads == 0 { available.min(MAX_THREADS) } else { requested_threads.min(available).min(MAX_THREADS) };

    let mut input = String::new();
    io::stdin().read_to_string(&mut input).unwrap_or_else(|error| fail(&format!("read stdin: {}", error)));
    match mode.as_str() {
        "align" => run_alignment(&input, threads),
        "disulfide" => run_disulfide(&input, threads),
        _ => fail("usage: protstruct-core <align|disulfide> [--threads N] < input.tsv"),
    }
}
