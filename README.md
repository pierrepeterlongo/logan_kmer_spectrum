README
# logan k-mer Spectrum 

This Rust program processes a FASTA file (including `.zst` compressed files) to compute a k-mer spectrum, where each k-mer's count is weighted by the abundance found in the FASTA headers, following the logan format. 

## Features
- Supports **both uncompressed and `.zst` compressed** FASTA files.
- Extracts **abundance** from FASTA headers in the format:
  ```
  >[accession]_[counter] ka:f:[abundance]
  ```
- Computes k-mer frequencies and their histogram. A k-mer count is provided by sum of its abundance in each sequence where it is found.
- Supports an **optional limit** to restrict output frequencies.

## Installation
Ensure you have Rust installed. If not, install it via [Rustup](https://rustup.rs/):
```sh
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Clone the repository and install dependencies:
```sh
git clone https://github.com/pierrepeterlongo/logan_kmer_spectrum.git
cd logan_kmer_spectrum
cargo build --release
```

## Usage
Run the tool with a FASTA file:
```sh
logan_kmer_spectrum --fasta_file input.fasta --k 31
```
For `.zst` compressed FASTA files:
```sh
logan_kmer_spectrum --fasta_file input.fasta.zst --k 31
```
To **limit** the maximum frequency displayed:
```sh
logan_kmer_spectrum --fasta_file input.fasta --k 31 --limit 100
```

## Dependencies
The following crates are used:
```toml
[dependencies]
bio = "0.41"
clap = { version = "4.4", features = ["derive"] }
regex = "1.10"
zstd = "0.13"
```

## Output
The program prints a frequency histogram:
```
K-mer Frequency    Count
1    20
2    35
3    50
...
```
If `--limit` is set, it restricts the output to that max frequency.

## License
MIT License.

## Author
Your Name

## Contributions
Pull requests are welcome! Open an issue for bug reports or feature requests.

# logan kmer spectrum

This Rust program processes a FASTA file (including `.zst` compressed files) to compute a k-mer spectrum, where each k-mer's count is weighted by the abundance found in the FASTA headers. The k-mer count is the sum of its abundances of sequences to which it belongs to. 

## Features
- Supports **both uncompressed and `.zst` compressed** FASTA files.
- Extracts **abundance** from FASTA headers in the format:
  ```
  >[accession]_[counter] ka:f:[abundance]
  ```
- Computes k-mer frequencies and their histogram.
- Supports an **optional limit** to restrict output frequencies.
- Optionally considers all k-mers as **canonical** (i.e., the lexicographically smallest representation between a k-mer and its reverse complement).

## Installation
Ensure you have Rust installed. If not, install it via [Rustup](https://rustup.rs/):
```sh
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Clone the repository and install dependencies:
```sh
git clone https://github.com/your-repo/kmer-spectrum-tool.git
cd kmer-spectrum-tool
cargo build --release
```

## Usage
Run the tool with a FASTA file:
```sh
cargo run --release -- --fasta_file input.fasta --k 31
```
For `.zst` compressed FASTA files:
```sh
cargo run --release -- --fasta_file input.fasta.zst --k 31
```
To **limit** the maximum frequency displayed:
```sh
cargo run --release -- --fasta_file input.fasta --k 31 --limit 100
```

## Output
The program prints a frequency histogram:
```
K-mer_Frequency    Count
1    20
2    35
3    50
5    10
...
```
If `--limit` is set, it restricts the output to that max frequency.

## License
AGP-L 3

## Author
Pierre Peterlongo

## Contributions
Pull requests are welcome! Open an issue for bug reports or feature requests.

