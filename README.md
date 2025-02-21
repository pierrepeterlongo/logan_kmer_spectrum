# logan kmer spectrum

This Rust program processes a FASTA file (including `.zst` compressed files) to compute a k-mer spectrum, where each k-mer's count is weighted by the abundance found in the FASTA headers, using the [logan format](https://github.com/IndexThePlanet/Logan). The k-mer count is the sum of its abundances in sequences to which it belongs to. 

## Features
- Supports **both uncompressed and `.zst` compressed** FASTA files.
- Extracts **abundance** from FASTA headers in the format:
  ```
  >[accession]_[counter] ka:f:[abundance]
  ```
- Computes k-mer frequencies and their histogram.
- Supports an **optional limit** to restrict output frequencies.
- Optionally considers all k-mers as **canonical** (i.e., the lexicographically smallest representation between a k-mer and its reverse complement).
- Not parallelized
- Not memory efficient


## Installation
Ensure you have Rust installed. If not, install it via [Rustup](https://rustup.rs/):
```sh
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Clone the repository and install dependencies:
```sh
git clone https://github.com/your-repo/logan_kmer_spectrum.git
cd logan_kmer_spectrum
cargo install --path .
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
To consider all k-mers as canonical:
```
logan_kmer_spectrum --fasta_file input.fasta --k 31 --canonical
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

