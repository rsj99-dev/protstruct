//! PSIPRED V4.1 单序列推理库，供命令行和批量检索核心复用。

mod blosum;

use blosum::BLOSUM62;

const P1_IPERGRP: usize = 41;
const P1_WINR: isize = 17;
const P1_CWIDTH: usize = 9;
const P1_CDEPTH: usize = 45;
const P1_NUM_IN: usize = 35 * P1_IPERGRP;
const P1_NUM_CONV: usize = P1_CDEPTH * (35 - P1_CWIDTH + 1);
const P1_NUM_HID: usize = 40;
const P1_NUM_OUT: usize = 3;
const P1_TOTAL: usize = P1_NUM_IN + P1_NUM_CONV + P1_NUM_HID + P1_NUM_OUT;
const P1_CONV_W: usize = P1_CWIDTH * P1_IPERGRP;

const P2_IPERGRP: usize = 4;
const P2_WINR: isize = 7;
const P2_NUM_IN: usize = 15 * P2_IPERGRP + 3 + 1;
const P2_NUM_HID: usize = 25;
const P2_NUM_OUT: usize = 3;
const P2_TOTAL: usize = P2_NUM_IN + P2_NUM_HID + P2_NUM_OUT;
const P2_OFFSET: usize = 15 * P2_IPERGRP;

const AA_ORDER: [u8; 20] = *b"ACDEFGHIKLMNPQRSTVWY";

struct Pass1Weights {
    conv: Vec<f32>,
    hidden: Vec<f32>,
    out: Vec<f32>,
    bias: Vec<f32>,
}

struct Pass2Weights {
    hidden: Vec<f32>,
    out: Vec<f32>,
    bias: Vec<f32>,
}

/// 一次加载权重后可重复预测，避免批量任务重复读取模型。
pub struct PsipredModel {
    pass1: [Pass1Weights; 3],
    pass2: Pass2Weights,
}

fn read_numbers(path: &str) -> Result<Vec<f32>, String> {
    let text = std::fs::read_to_string(path).map_err(|error| format!("read {}: {}", path, error))?;
    let mut values = Vec::with_capacity(text.len() / 4);
    for token in text.split_whitespace() {
        let value = token.parse::<f32>().map_err(|_| format!("bad weights file: {}", path))?;
        values.push(value);
    }
    Ok(values)
}

fn load_pass1(path: &str) -> Result<Pass1Weights, String> {
    let values = read_numbers(path)?;
    let conv_count = P1_NUM_CONV * P1_CONV_W;
    let hidden_count = P1_NUM_HID * P1_NUM_CONV;
    let out_count = P1_NUM_OUT * P1_NUM_HID;
    let bias_count = P1_NUM_CONV + P1_NUM_HID + P1_NUM_OUT;
    let expected = conv_count + hidden_count + out_count + bias_count;
    if values.len() != expected {
        return Err(format!("unexpected PSIPRED pass-1 weight count: {} (expected {})", values.len(), expected));
    }
    let mut offset = 0;
    let conv = values[offset..offset + conv_count].to_vec();
    offset += conv_count;
    let hidden = values[offset..offset + hidden_count].to_vec();
    offset += hidden_count;
    let out = values[offset..offset + out_count].to_vec();
    offset += out_count;
    let bias = values[offset..offset + bias_count].to_vec();
    Ok(Pass1Weights { conv, hidden, out, bias })
}

fn load_pass2(path: &str) -> Result<Pass2Weights, String> {
    let values = read_numbers(path)?;
    let hidden_count = P2_NUM_HID * P2_NUM_IN;
    let out_count = P2_NUM_OUT * P2_NUM_HID;
    let bias_count = P2_NUM_HID + P2_NUM_OUT;
    let expected = hidden_count + out_count + bias_count;
    if values.len() != expected {
        return Err(format!("unexpected PSIPRED pass-2 weight count: {} (expected {})", values.len(), expected));
    }
    let mut offset = 0;
    let hidden = values[offset..offset + hidden_count].to_vec();
    offset += hidden_count;
    let out = values[offset..offset + out_count].to_vec();
    offset += out_count;
    let bias = values[offset..offset + bias_count].to_vec();
    Ok(Pass2Weights { hidden, out, bias })
}

impl PsipredModel {
    pub fn load(data_dir: &str) -> Result<Self, String> {
        Ok(Self {
            pass1: [
                load_pass1(&format!("{}/weights.dat", data_dir))?,
                load_pass1(&format!("{}/weights.dat2", data_dir))?,
                load_pass1(&format!("{}/weights.dat3", data_dir))?,
            ],
            pass2: load_pass2(&format!("{}/weights_p2.dat", data_dir))?,
        })
    }

    pub fn predict(&self, sequence: &str) -> Result<String, String> {
        let sequence = sequence.trim().to_ascii_uppercase();
        if sequence.len() < 5 {
            return Err("PSIPRED requires a sequence of at least 5 residues".to_string());
        }
        if sequence.len() > 10000 {
            return Err("sequence is too long for PSIPRED".to_string());
        }
        let scores = pass1_scores(sequence.as_bytes(), &self.pass1);
        Ok(pass2_states(sequence.as_bytes(), &scores, &self.pass2))
    }
}

fn aa_num(ch: u8) -> usize {
    const MAP: [usize; 27] = [
        999, 0, 20, 4, 3, 6, 13, 7, 8, 9, 20, 11, 10, 12, 2, 20, 14, 5, 1, 15, 16, 20, 19, 17, 20, 18, 20,
    ];
    let upper = ch.to_ascii_uppercase();
    if upper.is_ascii_alphabetic() { MAP[(upper & 31) as usize] } else { 20 }
}

fn aa_num_mtx(ch: u8) -> usize {
    const MAP: [usize; 27] = [
        999, 0, 20, 4, 3, 6, 13, 7, 8, 9, 22, 11, 10, 12, 2, 22, 14, 5, 1, 15, 16, 22, 19, 17, 22, 18, 21,
    ];
    let upper = ch.to_ascii_uppercase();
    if upper.is_ascii_alphabetic() { MAP[(upper & 31) as usize] } else { 22 }
}

fn logistic(value: f32) -> f32 {
    1.0 / (1.0 + (-value).exp())
}

fn relu(value: f32) -> f32 {
    if value < 0.0 { 0.0 } else { value }
}

fn profile_from_sequence(sequence: &[u8]) -> Vec<[i32; 20]> {
    let mut profile = Vec::with_capacity(sequence.len());
    for &residue in sequence {
        let row = aa_num_mtx(residue);
        let mut values = [0i32; 20];
        for &code in AA_ORDER.iter() {
            let column = aa_num_mtx(code);
            values[column] = BLOSUM62[row][column] as i32 * 100;
        }
        profile.push(values);
    }
    profile
}

fn pass1_scores(sequence: &[u8], weights: &[Pass1Weights; 3]) -> Vec<[f32; 3]> {
    let length = sequence.len();
    let profile = profile_from_sequence(sequence);
    let mut averaged = vec![[0.0f32; 3]; length];
    let mut confidence_sum = vec![0.0f32; length];
    let mut activation = vec![0.0f32; P1_TOTAL];

    for weight_set in weights {
        for position in 0..length {
            activation.fill(0.0);
            for offset in -P1_WINR..=P1_WINR {
                let window_index = (offset + P1_WINR) as usize;
                let base = window_index * P1_IPERGRP;
                let absolute = position as isize + offset;
                if absolute >= 0 && (absolute as usize) < length {
                    let residue = absolute as usize;
                    for aa in 0..20 {
                        activation[base + aa] = profile[residue][aa] as f32 / 1000.0;
                    }
                    let identity = aa_num(sequence[residue]);
                    if identity < 20 {
                        activation[base + 20 + identity] = 1.0;
                    } else {
                        activation[base + 40] = 1.0;
                    }
                } else {
                    activation[base + 40] = 1.0;
                }
            }

            for unit in 0..P1_NUM_CONV {
                let start = (unit * P1_IPERGRP) / P1_CDEPTH;
                let mut value = weight_set.bias[unit];
                let weight_offset = unit * P1_CONV_W;
                for j in 0..P1_CONV_W {
                    value += activation[start + j] * weight_set.conv[weight_offset + j];
                }
                activation[P1_NUM_IN + unit] = relu(value);
            }
            for unit in 0..P1_NUM_HID {
                let mut value = weight_set.bias[P1_NUM_CONV + unit];
                let weight_offset = unit * P1_NUM_CONV;
                for previous in 0..P1_NUM_CONV {
                    value += activation[P1_NUM_IN + previous] * weight_set.hidden[weight_offset + previous];
                }
                activation[P1_NUM_IN + P1_NUM_CONV + unit] = relu(value);
            }
            for unit in 0..P1_NUM_OUT {
                let mut value = weight_set.bias[P1_NUM_CONV + P1_NUM_HID + unit];
                let weight_offset = unit * P1_NUM_HID;
                for previous in 0..P1_NUM_HID {
                    value += activation[P1_NUM_IN + P1_NUM_CONV + previous] * weight_set.out[weight_offset + previous];
                }
                activation[P1_TOTAL - P1_NUM_OUT + unit] = logistic(value);
            }
            let c = activation[P1_TOTAL - 3];
            let h = activation[P1_TOTAL - 2];
            let e = activation[P1_TOTAL - 1];
            let confidence = 2.0 * c.max(h).max(e) + c.min(h).min(e) - c - h - e;
            averaged[position][0] += confidence * c;
            averaged[position][1] += confidence * h;
            averaged[position][2] += confidence * e;
            confidence_sum[position] += confidence;
        }
    }

    for position in 0..length {
        let denominator = confidence_sum[position];
        if denominator != 0.0 {
            averaged[position][0] /= denominator;
            averaged[position][1] /= denominator;
            averaged[position][2] /= denominator;
        }
    }
    averaged
}

fn round3(value: f32) -> f32 {
    (value * 1000.0).round() / 1000.0
}

fn pass2_states(sequence: &[u8], scores: &[[f32; 3]], weights: &Pass2Weights) -> String {
    let length = sequence.len();
    let mut profile = vec![[0.0f32; 3]; length];
    for position in 0..length {
        profile[position][0] = round3(scores[position][0]);
        profile[position][1] = round3(scores[position][1]);
        profile[position][2] = round3(scores[position][2]);
    }

    let average_c = profile.iter().map(|value| value[0]).sum::<f32>() / length as f32;
    let average_h = profile.iter().map(|value| value[1]).sum::<f32>() / length as f32;
    let average_e = profile.iter().map(|value| value[2]).sum::<f32>() / length as f32;
    let length_log = (length as f32).ln();

    let mut activation = vec![0.0f32; P2_TOTAL];
    let mut states = vec![b'C'; length];
    let mut score_c = vec![0.0f32; length];
    let mut score_h = vec![0.0f32; length];
    let mut score_e = vec![0.0f32; length];

    for position in 0..length {
        activation.fill(0.0);
        activation[P2_OFFSET] = average_c;
        activation[P2_OFFSET + 1] = average_h;
        activation[P2_OFFSET + 2] = average_e;
        activation[P2_OFFSET + 3] = length_log;
        for offset in -P2_WINR..=P2_WINR {
            let window_index = (offset + P2_WINR) as usize;
            let base = window_index * P2_IPERGRP;
            let absolute = position as isize + offset;
            if absolute >= 0 && (absolute as usize) < length {
                let residue = absolute as usize;
                activation[base] = profile[residue][0];
                activation[base + 1] = profile[residue][1];
                activation[base + 2] = profile[residue][2];
            } else {
                activation[base + 3] = 1.0;
            }
        }
        for unit in 0..P2_NUM_HID {
            let mut value = weights.bias[unit];
            let weight_offset = unit * P2_NUM_IN;
            for input in 0..P2_NUM_IN {
                value += activation[input] * weights.hidden[weight_offset + input];
            }
            activation[P2_NUM_IN + unit] = relu(value);
        }
        for unit in 0..P2_NUM_OUT {
            let mut value = weights.bias[P2_NUM_HID + unit];
            let weight_offset = unit * P2_NUM_HID;
            for previous in 0..P2_NUM_HID {
                value += activation[P2_NUM_IN + previous] * weights.out[weight_offset + previous];
            }
            activation[P2_TOTAL - P2_NUM_OUT + unit] = logistic(value);
        }
        let c = activation[P2_TOTAL - 3];
        let h = activation[P2_TOTAL - 2];
        let e = activation[P2_TOTAL - 1];
        states[position] = if c > h && c > e {
            b'C'
        } else if h > c && h > e {
            b'H'
        } else {
            b'E'
        };
        score_c[position] = c;
        score_h[position] = h;
        score_e[position] = e;
    }

    let confidence: Vec<f32> = (0..length)
        .map(|position| {
            let c = score_c[position];
            let h = score_h[position];
            let e = score_e[position];
            2.0 * c.max(h).max(e) + c.min(h).min(e) - c - h - e
        })
        .collect();

    // PSIPRED pass-2 singleton helix/strand filter.
    for position in 1..length.saturating_sub(1) {
        if states[position] != b'C'
            && states[position - 1] == states[position + 1]
            && confidence[position] < 0.5 * (confidence[position - 1] + confidence[position + 1])
        {
            states[position] = states[position - 1];
        }
    }
    // PSIPRED pass-2 terminal C filter.
    for position in 1..length.saturating_sub(1) {
        if states[position - 1] == b'C' && states[position] != states[position + 1] {
            states[position] = b'C';
        }
        if states[position + 1] == b'C' && states[position] != states[position - 1] {
            states[position] = b'C';
        }
    }

    String::from_utf8(states).expect("internal state encoding error")
}
