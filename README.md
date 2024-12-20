# Comparative Genomics with BLAST

This repository contains a custom implementation of the **Basic Local Alignment Search Tool (BLAST)**, tailored for identifying homologous genomic regions between species using computationally efficient methods. This work focuses on adapting the BLAST algorithm to run effectively on storage-limited systems while maintaining accuracy and biological relevance.

## Overview
Comparative genomics is a field of research that compares the genomes of different species to understand evolutionary relationships, conserved genetic elements, and species-specific adaptations. Our implementation of BLAST is designed to identify homologous regions between genomic sequences, focusing on:

- Human and chimpanzee genomes.
- Selected chromosomes: 7 and 20.
- Reduced computational complexity suitable for low-end systems.

## Features
- **Seed and Extend Algorithm**: Uses k-mers (short sequences of fixed length) to identify regions of local similarity.
- **Custom Modifications**:
  - Selects k-mers with moderate matches to reduce noise and improve efficiency.
  - Extends alignments using a scoring system based on match and mismatch penalties.
  - Threshold-based termination of extension to avoid excessive computation.
- **Global Alignment Validation**: Confirms local alignments using a custom global alignment algorithm based on Fast Fourier Transform (FFT).
- **Visualization Tools**:
  - Plot alignment scores and frequencies.
  - Highlight homologous regions graphically.

## Dependencies
- Python 3.8+
- `numpy`
- `scipy`
- `matplotlib`
- `Biopython`

## File Structure
- **`blast.py`**: Core BLAST algorithm implementation.
- **`global_alignment.py`**: Implements FFT-based global alignment to validate homologous regions.
- **`visualization.py`**: Functions to visualize alignment results.
- **`utils.py`**: Utility functions for handling FASTA files and sequence preprocessing.

## Getting Started

### 1. Installation
Ensure you have Python installed. Install the required dependencies:
```bash
pip install numpy scipy matplotlib biopython
```

### 2. Usage

#### Input
Provide DNA sequences as FASTA files. Use the `fasta_to_str` utility to load sequences:
```python
from utils import fasta_to_str
seq1 = fasta_to_str("path/to/seq1.fasta")
seq2 = fasta_to_str("path/to/seq2.fasta")
```

#### BLAST Algorithm
Initialize the BLAST class and perform alignment:
```python
from blast import BLAST
blast = BLAST(seq1, seq2)
blast.seed(kmers=10, kmer_length=10)
blast.extend_seeds(perc_homology=0.9, min_match_len=20000, threshold=1000)
blast.print_alignment_report()
```

#### Visualization
Plot alignment results:
```python
from visualization import plot_alignment, plot_frequency
plot_alignment(seq1, blast.matches, "Alignment Plot")
plot_frequency(blast.freq, "Frequency Plot")
```

## Results
Our modified BLAST algorithm successfully identifies conserved regions between human and chimpanzee genomes, with:
- High accuracy for local alignments.
- Computational efficiency on chromosomes 7 and 20.

## Limitations
- Less effective for species with significant evolutionary divergence (e.g., human vs. mouse).
- No support for structural variations such as large insertions or deletions.

## Future Work
- Implementing gap penalties for handling structural variations.
- Enhancing scoring schemes for better alignment accuracy.
- Expanding support to more chromosomes and distant species comparisons.

## Citation
If you use this implementation in your work, please cite this repository

## License
This project is licensed under the MIT License.

---
For questions or contributions, please open an issue or submit a pull request.

