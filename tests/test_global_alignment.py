import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

import pytest
import numpy as np

from global_alignment import GlobalAlignment

def test_global_alignment_initialization():
    # Test initialization with valid sequences
    seq1 = "ACGTACGT"
    seq2 = "TACGTACG"
    alignment = GlobalAlignment(seq1, seq2)
    assert alignment.seq1 == seq1
    assert alignment.seq2 == seq2
    assert len(alignment.matches) == len(seq1) + len(seq2) - 1

    # Test initialization with empty sequences
    with pytest.raises(ValueError):
        GlobalAlignment("", "ACGT")
    with pytest.raises(ValueError):
        GlobalAlignment("ACGT", "")

def test_align_strands():
    seq1 = "ACGTACGT"
    seq2 = "TACGTACG"
    alignment = GlobalAlignment(seq1, seq2)

    matches = alignment.align_strands()
    assert isinstance(matches, np.ndarray)
    assert matches.shape[0] == len(seq1) + len(seq2) - 1

def test_matches_to_freq():
    seq1 = "ACGTACGT"
    seq2 = "TACGTACG"
    alignment = GlobalAlignment(seq1, seq2)

    freq = alignment.matches_to_freq()
    assert isinstance(freq, np.ndarray)
    assert freq.shape[0] == len(seq1) + len(seq2) - 1
    assert np.all(freq >= 0) and np.all(freq <= 1)

def test_max_alignment():
    seq1 = "ACGTACGT"
    seq2 = "TACGTACG"
    alignment = GlobalAlignment(seq1, seq2)

    max_index = alignment.max_alignment
    assert isinstance(max_index, int)
    assert 0 <= max_index < len(alignment.matches)

def test_print_report(capsys):
    seq1 = "ACGTACGT"
    seq2 = "TACGTACG"
    alignment = GlobalAlignment(seq1, seq2)

    alignment.print_report()
    captured = capsys.readouterr()
    assert "seq1 length" in captured.out
    assert "seq2 length" in captured.out
    assert "Most bases aligned" in captured.out
    assert "Index to align" in captured.out
