//! PSIPRED V4.1 单序列与批量 FASTA 推理命令行。

use std::env;
use std::fs;
use std::process;

use protstruct_psipred::PsipredModel;

fn fail(message: &str) -> ! {
    eprintln!("{}", message);
    process::exit(2);
}

fn parse_fasta(text: &str) -> Vec<(String, String)> {
    let mut records = Vec::new();
    let mut identifier: Option<String> = None;
    let mut sequence = String::new();
    for raw_line in text.lines() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with(';') || line.starts_with('#') {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            if let Some(previous) = identifier.take() {
                records.push((previous, std::mem::take(&mut sequence)));
            }
            let header = header.trim();
            if header.is_empty() {
                fail("empty FASTA header");
            }
            identifier = Some(header.split_whitespace().next().unwrap_or("sequence").to_string());
        } else if identifier.is_some() {
            for ch in line.chars() {
                if !ch.is_whitespace() && ch != '-' && ch != '.' && ch != '*' {
                    sequence.push(ch.to_ascii_uppercase());
                }
            }
        } else {
            fail("sequence data found before the first '>' header");
        }
    }
    if let Some(previous) = identifier {
        records.push((previous, sequence));
    }
    if records.is_empty() {
        fail("no sequences found");
    }
    records
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let mut data_dir = String::new();
    let mut sequence = String::new();
    let mut input = String::new();
    let mut output = String::new();
    let mut index = 1;
    while index < args.len() {
        match args[index].as_str() {
            "--data" if index + 1 < args.len() => {
                data_dir = args[index + 1].clone();
                index += 2;
            }
            "--sequence" if index + 1 < args.len() => {
                sequence = args[index + 1].clone();
                index += 2;
            }
            "--input" if index + 1 < args.len() => {
                input = args[index + 1].clone();
                index += 2;
            }
            "--output" if index + 1 < args.len() => {
                output = args[index + 1].clone();
                index += 2;
            }
            "--help" | "-h" => {
                println!("usage: protstruct-psipred --data <weights-dir> (--sequence <SEQ> | --input <file.fa>) [--output <file>]");
                return;
            }
            value if sequence.is_empty() && input.is_empty() => {
                sequence = value.to_string();
                index += 1;
            }
            _ => fail("unknown or incomplete argument"),
        }
    }
    if data_dir.is_empty() || (sequence.is_empty() && input.is_empty()) || (!sequence.is_empty() && !input.is_empty()) {
        fail("usage: protstruct-psipred --data <weights-dir> (--sequence <SEQ> | --input <file.fa>) [--output <file>]");
    }

    let model = match PsipredModel::load(&data_dir) {
        Ok(model) => model,
        Err(error) => fail(&error),
    };
    let result = if sequence.is_empty() {
        let text = fs::read_to_string(&input).unwrap_or_else(|error| fail(&format!("read {}: {}", input, error)));
        let records = parse_fasta(text.trim_start_matches('\u{feff}'));
        let mut rendered = String::new();
        for (identifier, residues) in records {
            match model.predict(&residues) {
                Ok(states) => {
                    rendered.push_str(&format!(">{}\n{}\n", identifier, states));
                }
                Err(error) => fail(&format!("{}: {}", identifier, error)),
            }
        }
        rendered
    } else {
        match model.predict(&sequence) {
            Ok(states) => states,
            Err(error) => fail(&error),
        }
    };

    if output.is_empty() {
        print!("{}", result);
        if !result.ends_with('\n') {
            println!();
        }
    } else if let Err(error) = fs::write(&output, result) {
        fail(&format!("write {}: {}", output, error));
    }
}
