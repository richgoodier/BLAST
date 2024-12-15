import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

import pytest
from blast import BLAST

def test_blast_initialization():
    """Test if the BLAST class initializes correctly."""
    seq1 = "ACTGACTGACTG"
    seq2 = "ACTGACTGACCC"
    blast = BLAST(seq1, seq2)
    assert blast.seq1 == seq1
    assert blast.seq2 == seq2
    assert blast.possible_kmers == []
    assert blast.seq1_index == {}
    assert blast.seq2_index == {}
    assert blast.alignment_info_list == []

def test_seed_generation():
    """Test the seed generation functionality."""
    seq1 = "ACTGACTGACTG"
    seq2 = "ACTGACTGACCC"
    blast = BLAST(seq1, seq2)
    blast.seed(kmers=5, kmer_length=4, verbose=False)
    assert len(blast.possible_kmers) <= 5
    for kmer in blast.possible_kmers:
        assert len(kmer) == 4
        assert kmer in blast.seq1_index
        assert kmer in blast.seq2_index

def test_seed_invalid_kmer_length():
    """Test seed generation with invalid k-mer length."""
    seq1 = "ACTG"
    seq2 = "ACTG"
    blast = BLAST(seq1, seq2)
    blast.seed(kmer_length=10, verbose=False)
    assert blast.possible_kmers == []

def test_extend_seeds():
    """Test extending the seeds functionality."""
    seq1 = "ACTGACTGACTGACTGACTG"
    seq2 = "ACTGACTGACTGCCCCACTG"
    blast = BLAST(seq1, seq2)
    blast.seed(kmers=3, kmer_length=4, verbose=False)
    blast.extend_seeds(perc_homology=0.8, min_match_len=8, threshold=5)
    assert len(blast.alignment_info_list) > 0
    for alignment in blast.alignment_info_list:
        assert alignment['match_len'] >= 8

def test_homology_calculation():
    """Test the homology calculation method."""
    seq1 = "ACTGACTG"
    seq2 = "ACTGACCG"
    blast = BLAST(seq1, seq2)
    matches, percent_homology = blast.homology(seq1, seq2)
    assert matches == 7
    assert percent_homology == pytest.approx(0.875)

def test_export_homologous_sequences():
    """Test exporting homologous sequences."""
    seq1 = "ACTGACTGACTGACTGACTG"
    seq2 = "ACTGACTGACTGCCCCACTG"
    blast = BLAST(seq1, seq2)
    blast.seed(kmers=3, kmer_length=4, verbose=False)
    blast.extend_seeds(perc_homology=0.8, min_match_len=8, threshold=5)

    if blast.alignment_info_list:
        kmer = blast.alignment_info_list[0]['kmer']
        seq1_start = blast.alignment_info_list[0]['seq1_start']
        seq2_start = blast.alignment_info_list[0]['seq2_start']
        homologous_sequences = blast.export_homologous_sequences(kmer, seq1_start, seq2_start)
        assert homologous_sequences is not None
    else:
        pytest.skip("No alignments found for testing export.")

def test_print_alignment_report():
    """Test printing alignment reports."""
    seq1 = "ACTGACTGACTGACTGACTG"
    seq2 = "ACTGACTGACTGCCCCACTG"
    blast = BLAST(seq1, seq2)
    blast.seed(kmers=3, kmer_length=4, verbose=False)
    blast.extend_seeds(perc_homology=0.8, min_match_len=8, threshold=5)
    blast.print_alignment_report()
    assert True  # Just ensure no exceptions during the process
