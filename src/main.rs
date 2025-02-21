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
    /// k-mer size (maximum 32 for 64-bit representation)
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

/// Convert a nucleotide to its 2-bit representation
const fn nucleotide_to_bits(n: u8) -> Option<u64> {
    match n {
        b'A' | b'a' => Some(0b00),
        b'C' | b'c' => Some(0b01),
        b'G' | b'g' => Some(0b10),
        b'T' | b't' => Some(0b11),
        _ => None,
    }
}

/// Convert a nucleotide to its complement's 2-bit representation
const fn complement_to_bits(n: u8) -> Option<u64> {
    match n {
        b'A' | b'a' => Some(0b11), // T
        b'C' | b'c' => Some(0b10), // G
        b'G' | b'g' => Some(0b01), // C
        b'T' | b't' => Some(0b00), // A
        _ => None,
    }
}

/// Encodes a DNA sequence into a 64-bit integer
fn encode_kmer(seq: &[u8], k: usize) -> Option<u64> {
    if k > 32 || seq.len() < k {
        return None;
    }
    
    let mut encoded: u64 = 0;
    for i in 0..k {
        if let Some(bits) = nucleotide_to_bits(seq[i]) {
            encoded = (encoded << 2) | bits;
        } else {
            return None; // Invalid nucleotide
        }
    }
    Some(encoded)
}

/// Encodes the reverse complement of a DNA sequence
fn encode_reverse_complement(seq: &[u8], k: usize) -> Option<u64> {
    if k > 32 || seq.len() < k {
        return None;
    }
    
    let mut encoded: u64 = 0;
    for i in (0..k).rev() {
        if let Some(bits) = complement_to_bits(seq[i]) {
            encoded = (encoded << 2) | bits;
        } else {
            return None; // Invalid nucleotide
        }
    }
    Some(encoded)
}

/// Generates k-mers as bit-encoded u64 values, considering canonical representation if required
fn generate_encoded_kmers(seq: &[u8], k: usize, canonical: bool) -> Vec<u64> {
    if seq.len() < k {
        return vec![];
    }

    let mut kmers: Vec<u64> = Vec::new();
    for i in 0..=seq.len() - k {
        let kmer_slice = &seq[i..i + k];
        
        if let Some(encoded) = encode_kmer(kmer_slice, k) {
            if canonical {
                if let Some(rev_comp) = encode_reverse_complement(kmer_slice, k) {
                    kmers.push(std::cmp::min(encoded, rev_comp));
                }
            } else {
                kmers.push(encoded);
            }
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
    
    if args.k > 32 {
        eprintln!("Error: k-mer size cannot exceed 32 for the 64-bit representation");
        std::process::exit(1);
    }
    
    let fasta_path = Path::new(&args.fasta_file);
    let fasta_reader = open_fasta_file(fasta_path)?;
    let reader = fasta::Reader::new(fasta_reader);
    
    let mut kmer_counts: HashMap<u64, u64> = HashMap::new();

    for result in reader.records() {
        let record = result?;
        let header = format!(
            "{} {}", 
            record.id(), 
            record.desc().unwrap_or("")
        );
        let sequence = record.seq();
        
        if let Some(abundance) = extract_abundance(&header) {
            let kmers = generate_encoded_kmers(sequence, args.k, args.canonical);
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