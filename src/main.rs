use bio::io::fasta;
use clap::Parser;
use regex::Regex;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;
use zstd::Decoder;

/// Command-line arguments
#[derive(Parser)]
struct Args {
    /// Input FASTA file (supports `.zst` compressed files)
    fasta_file: String,
    /// k-mer size
    k: usize,
    /// Optional: Maximum frequency to display
    #[arg(short, long)]
    limit: Option<u64>,
    /// Optional: Consider all k-mers as canonical
    #[arg(long)]
    canonical: bool,
}

/// Extracts abundance from FASTA header following `[accession]_[counter] ka:f:[abundance]`
fn extract_abundance(header: &str) -> Option<u32> {
    let re = Regex::new(r"ka:f:(\d+)").unwrap();
    re.captures(header).and_then(|caps| caps.get(1)).map(|m| m.as_str().parse().ok()).flatten()
}

/// Computes the reverse complement of a DNA sequence
fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => 'N'
        })
        .collect()
}

/// Generates k-mers, considering canonical representation if required
fn generate_kmers(seq: &str, k: usize, canonical: bool) -> Vec<String> {
    if seq.len() < k {
        return vec![];
    }

    let mut kmers: Vec<String> = Vec::new();
    for i in 0..=seq.len() - k {
        let kmer = seq[i..i + k].to_string();
        if canonical {
            let rev_kmer = reverse_complement(&kmer);
            kmers.push(std::cmp::min(kmer, rev_kmer));
        } else {
            kmers.push(kmer);
        }
    }
    kmers
}

/// Opens a FASTA file, supporting both regular and `.zst` compressed formats
fn open_fasta_file(file_path: &Path) -> std::io::Result<Box<dyn Read>> {
    let file = File::open(file_path)?;
    if file_path.extension().map_or(false, |ext| ext == "zst") {
        let decoder = Decoder::new(file)?;
        Ok(Box::new(BufReader::new(decoder)))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn main() -> std::io::Result<()> {
    let args = Args::parse();
    let fasta_path = Path::new(&args.fasta_file);
    let fasta_reader = open_fasta_file(fasta_path)?;
    let reader = fasta::Reader::new(fasta_reader);
    
    let mut kmer_counts: HashMap<String, u64> = HashMap::new();

    for result in reader.records() {
        let record = result?;
        let header = format!(
            "{} {}", 
            record.id(), 
            record.desc().unwrap_or("")
        );
        let sequence = String::from_utf8_lossy(record.seq());
        
        if let Some(abundance) = extract_abundance(&header) {
            let kmers = generate_kmers(&sequence, args.k, args.canonical);
            for kmer in kmers {
                *kmer_counts.entry(kmer).or_insert(0) += abundance as u64;
            }
        }
    }

    // Compute histogram
    let mut histogram: HashMap<u64, u64> = HashMap::new();
    for &count in kmer_counts.values() {
        *histogram.entry(count).or_insert(0) += 1;
    }

    // Sort the histogram for printing
    let mut sorted_histogram: Vec<_> = histogram.iter().collect();
    sorted_histogram.sort();

    // Print header
    println!("K-mer Frequency\tCount");

    // Print histogram, applying the optional `--limit`
    for (freq, count) in sorted_histogram {
        if let Some(limit) = args.limit {
            if *freq > limit {
                break;
            }
        }
        println!("{}\t{}", freq, count);
    }

    Ok(())
}
