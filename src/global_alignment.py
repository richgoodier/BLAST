import numpy as np
from scipy.signal import fftconvolve

class GlobalAlignment:
    """
    Class to perform global alignment of two DNA sequences using FFT-based convolution.

    Attributes:
    -----------
    seq1 : str
        The first DNA sequence.
    seq2 : str
        The second DNA sequence.
    matches : np.ndarray
        An array of match scores for each alignment shift.
    max_alignment : int
        The index of the highest match score.
    freq : np.ndarray
        The normalized frequency of matches across alignment shifts.

    Methods:
    --------
    align_strands():
        Computes match scores for all alignment shifts using FFT-based convolution.
    matches_to_freq():
        Converts raw match scores to normalized frequencies based on overlapping bases.
    print_report():
        Prints a summary report of the alignment analysis.
    """

    def __init__(self, seq1: str, seq2: str) -> None:
        """
        Initialize the GlobalAlignment class with two DNA sequences.

        Parameters:
        -----------
        seq1 : str
            The first DNA sequence.
        seq2 : str
            The second DNA sequence.
        """
        self.seq1 = seq1
        self.seq2 = seq2
        if len(seq1) == 0 or len(seq2) == 0:
            raise ValueError("Ensure both sequences are non-empty.")
        self.matches = self.align_strands()
        self.max_alignment = np.argmax(self.matches)
        self.freq = self.matches_to_freq()

    def align_strands(self) -> np.ndarray:
        """
        Align two DNA sequences using FFT convolution and return match scores.

        Returns:
        --------
        np.ndarray
            Array of integers representing match scores for each alignment shift.
        """
        seq1_arr = np.array(list(self.seq1))
        seq2_arr = np.array(list(self.seq2))

        # Identify unique characters in each sequence
        seq1_chars = set(seq1_arr)
        seq2_chars = set(seq2_arr)
        nucleotides = list(seq1_chars.union(seq2_chars))

        len0 = len(seq1_arr)
        len1 = len(seq2_arr)
        total_shifts = len0 + len1 - 1

        matches = np.zeros(total_shifts, dtype='uint32')

        # Compute cross-correlations for each nucleotide using FFT convolution
        for nucleotide in nucleotides:
            seq0_binary = (seq1_arr == nucleotide).astype(int)
            seq1_binary = (seq2_arr == nucleotide).astype(int)

            # Reverse seq1_binary because fftconvolve performs convolution, not correlation
            seq1_binary_reversed = seq1_binary[::-1]

            # Compute convolution using FFT
            conv_result = fftconvolve(seq0_binary, seq1_binary_reversed, mode='full')

            # Since the result may have small imaginary parts due to numerical errors, take the real part
            conv_result = np.round(np.real(conv_result)).astype('uint32')

            matches += conv_result

        return matches

    def matches_to_freq(self) -> np.ndarray:
        """
        Normalize match scores to frequencies based on overlapping base counts.

        Returns:
        --------
        np.ndarray
            Array of normalized match frequencies.
        """
        seq1_len = len(self.seq1)
        seq2_len = len(self.seq2)
        overlap = min(seq1_len, seq2_len)
        misc = abs(seq1_len - seq2_len)

        front = np.arange(1, overlap)
        middle = np.full(misc + 1, overlap)
        back = front[::-1]
        
        overlaps = np.concatenate([front, middle, back])

        return self.matches / overlaps

    def print_report(self) -> None:
        """
        Print a summary report of the alignment analysis.
        """
        if len(self.matches) > 0:
            print(f"seq1 length: {len(self.seq1)}  seq2 length: {len(self.seq2)}")
            print("Most bases aligned:", self.matches[self.max_alignment])
            print("Index to align:", self.max_alignment - len(self.seq1))
