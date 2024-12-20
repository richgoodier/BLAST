import numpy as np
import itertools


class BLAST:
    """
    A class to perform a simplified version of the Basic Local Alignment Search Tool (BLAST) algorithm.

    This implementation identifies regions of similarity between two sequences using k-mer generation and a seed-and-extend approach.

    Attributes:
    -----------
    seq1 : str
        The first sequence to be compared.
    seq2 : str
        The second sequence to be compared.
    possible_kmers : list
        A list of selected k-mers from seq1.
    seq1_index : dict
        Maps k-mers to their positions in seq1.
    seq2_index : dict
        Maps k-mers to their positions in seq2.
    alignment_info_list : list
        Stores information about extended alignments between seq1 and seq2.
    
    Methods:
    --------
    seed(kmers=10, kmer_length=10, verbose=True):
        Generates k-mers from seq1 and finds their matches in seq2.
    
    extend_seeds(perc_homology=0.9, min_match_len=20000, threshold=1000):
        Extends alignments for k-mers between two sequences using a simplified BLAST approach.
    
    homology(splice1, splice2):
        Calculates the number of matches and the percentage of homology between two aligned sequences.
    
    print_alignment_report():
        Prints a detailed report of the alignments, sorted from longest to shortest.
    
    export_homologous_sequences(kmer, seq1_start, seq2_start):
        Exports the homologous sequences associated with a specific k-mer and start positions.
    """

    def __init__(self, seq1: str, seq2: str) -> None:
        """
        Initializes the BLAST class with two sequences to compare.

        Parameters:
        -----------
        seq1 : str
            The first sequence for comparison.
        seq2 : str
            The second sequence for comparison.
        """
        self.seq1 = seq1
        self.seq2 = seq2
        self.possible_kmers = []
        self.seq1_index = {}
        self.seq2_index = {}
        self.alignment_info_list = []

    def seed(self, kmers: int = 10, kmer_length: int =10, verbose: bool = True) -> None:
        """
        Generate k-mers from seq1 and find their matches in seq2.

        Parameters:
        -----------
        kmers : int, optional
            The number of k-mers to randomly generate (default is 10).
        kmer_length : int, optional
            The length of each k-mer (default is 10).
            If kmer_length is longer than either seq1 or seq2, the function will exit.
        verbose : bool, optional
            If True, prints information about the selected k-mers and their frequencies (default is True).

        Notes:
        ------
        - Only k-mers with more than 80% valid nucleotide characters (G, A, T, C) are considered.
          This allows for masking (N, g, a, t, c)
        - Avoids overly common or rare k-mers to improve alignment efficiency.
        """
        
        seq1_len = len(self.seq1)
        seq2_len = len(self.seq2)
        if kmer_length > seq1_len or kmer_length > seq2_len:
            print("Error: kmer_length longer than seq1 and/or seq2")
            return
        if kmer_length > seq1_len * 0.1 or kmer_length > seq2_len * 0.1:
            print("Warning: kmer_length longer than 10% of seq1 and/or seq2.")
        
        if self.possible_kmers:
            query = input("Possible kmers already exist.  Overwrite and continue? (y/n): ").strip().lower()
            if query not in ["y", "yes"]:
                print("Operation canceled.")
                return
            else:
                while True:
                    query = input("Extend or overwrite? (e/o): ").strip().lower()
                    if query == 'o':
                        self.possible_kmers = []
                        self.seq1_index = {}
                        self.seq2_index = {}
                        total_kmers = kmers
                        break
                    elif query == 'e':
                        total_kmers = kmers + len(self.possible_kmers)
                        break
        else:
            total_kmers = kmers
        
        rejects = 0
        
        while len(self.possible_kmers) < total_kmers:
            start = np.random.randint(seq1_len - kmer_length)
            kmer = self.seq1[start:start+kmer_length]
            valid_base_perc = sum(kmer.count(base) for base in "GATC") / kmer_length

            if valid_base_perc > 0.8:  # Eliminates any kmer that contains less than 80% GATC characters.
                match_idx = []
                for i in range(seq2_len - kmer_length):  # search for the kmer in seq2
                    if self.seq2[i:i+kmer_length] == kmer:
                        match_idx.append(i)

                if 0 < len(match_idx) < 50:  # Eliminates kmers that are too common
                    self.possible_kmers.append(kmer)
                    self.seq2_index[kmer] = match_idx
                    if verbose:
                        print(f"{kmer} ---- {len(match_idx)} ---- ", end='')

                    match_idx = []
                    for i in range(seq1_len - kmer_length):
                        if self.seq1[i:i+kmer_length] == kmer:
                            match_idx.append(i)
                    self.seq1_index[kmer] = match_idx
                    if verbose:
                        print(len(match_idx))
                    else:
                        print(len(self.possible_kmers), end='-', flush=True)
                else:
                    rejects += 1
            else:
                rejects += 1
        
        print(f"{len(self.possible_kmers)} kmers stored.")
        print(f"{rejects} total rejects.")
    
    def extend_seeds(self, perc_homology: float = 0.9, min_match_len: int = 20_000, threshold: int = 1_000):
        """
        Extend alignments for k-mers between two sequences using a simplified BLAST approach.

        Parameters:
        -----------
        perc_homology : float, optional
            The desired percentage of homology (0.0 to 1.0) for alignment (default is 0.9).
            True homology also drops as threshold value increases.
        min_match_len : int, optional
            The minimum match length required to save an alignment (default is 20,000).
        threshold : int, optional
            The score threshold for terminating alignment extension (default is 1,000).

        Notes:
        ------
        - Alignments are extended in both directions from k-mer positions.
        - The algorithm stops extending when the threshold score is reached.
        """
        
        if not 0 <= perc_homology <= 1:
            print("Error: perc_homology must be between 0 and 1, inclusive.")
            return
        
        if self.alignment_info_list:
            query = input("Alignment info already exist.  Overwrite and continue? (y/n): ").strip().lower()
            if query not in ["y", "yes"]:
                print("Operation canceled.")
                return
            else:
                self.alignment_info_list = []
        
        for kmer in self.possible_kmers:
            for kmer_position_seq1, kmer_position_seq2 in itertools.product(self.seq1_index[kmer], self.seq2_index[kmer]):
                alignment_info = {}
                alignment_info["kmer"] = kmer
                match_len = len(kmer)

                # Extend right
                remaining_threshold = threshold
                i = kmer_position_seq1 + len(kmer)
                j = kmer_position_seq2 + len(kmer)
                while i < len(self.seq1) and j < len(self.seq2) and remaining_threshold > 0:
                    if self.seq1[i] == self.seq2[j]:
                        remaining_threshold += (1 - perc_homology)
                    else:
                        remaining_threshold -= 1
                    i += 1
                    j += 1
                    match_len += 1
                alignment_info['seq1_end'] = i
                alignment_info['seq2_end'] = j

                # Extend left
                remaining_threshold = threshold
                i = kmer_position_seq1 - 1
                j = kmer_position_seq2 - 1
                while i >= 0 and j >= 0 and remaining_threshold > 0:
                    if self.seq1[i] == self.seq2[j]:
                        remaining_threshold += (1 - perc_homology)
                    else:
                        remaining_threshold -= 1
                    i -= 1
                    j -= 1
                    match_len += 1
                alignment_info['seq1_start'] = i
                alignment_info['seq2_start'] = j

                alignment_info['match_len'] = match_len

                if alignment_info['match_len'] > min_match_len:
                    if alignment_info not in self.alignment_info_list:
                        self.alignment_info_list.append(alignment_info)
                    print(f"kmer: {alignment_info['kmer']}  match length: {alignment_info['match_len']}")
                    print(f"seq1[{alignment_info['seq1_start']}:{alignment_info['seq1_end']}]  seq2[{alignment_info['seq2_start']}:{alignment_info['seq2_end']}]")
        
        self.alignment_info_list = sorted(self.alignment_info_list, key=lambda x: x['match_len'], reverse=True)
        print(f"Recorded a total of {len(self.alignment_info_list)} homologous regions greater than {min_match_len} bases long.")
            
    def homology(self, splice1 : str, splice2 : str) -> tuple:
        """
        Calculate the number of matches and the percentage of homology between two aligned sequences.

        Parameters:
        -----------
        splice1 : str
            The first sequence segment.
        splice2 : str
            The second sequence segment.

        Returns:
        --------
        tuple
            A tuple containing the number of matches and the percentage of homology.
        
        Example:
        --------
            splice1 = "ACTGACTG"
            splice2 = "ACTGACCG"
            output: (7, 0.875)  # 7 matches, 87.5% homology
        """
        if len(splice1) != len(splice2):
            print("Error: Sequence lengths must be identical.")
            return 0, 0.0
        
        length_seq = len(splice1)
        if length_seq == 0:
            return 0, 0.0  # Handle empty sequences.
        
        num_matches = 0
        for i in range(length_seq):
            if splice1[i] == splice2[i]:
                num_matches += 1
        
        return num_matches, num_matches / length_seq
        
    def print_alignment_report(self, top_candidates = np.inf) -> None:
        """
        Print a detailed report of the alignments.

        The report includes the k-mer, match length, alignment positions, 
        and homology percentage for each alignment.
        """
        
        print("From longest to shortest:")
        i = 1
        for alignment_info in self.alignment_info_list:
            print(f"kmer: {alignment_info['kmer']}  match length: {alignment_info['match_len']}")
            print(f"seq1[{alignment_info['seq1_start']}:{alignment_info['seq1_end']}]  seq2[{alignment_info['seq2_start']}:{alignment_info['seq2_end']}]")
            splice1 = self.seq1[alignment_info['seq1_start']:alignment_info['seq1_end']]
            splice2 = self.seq2[alignment_info['seq2_start']:alignment_info['seq2_end']]
            num_matches, percent_homology = self.homology(splice1, splice2)
            print(f"A total of {num_matches} of {len(splice1)} bases in alignment.")
            print(f"Fraction aligned: {percent_homology}")
            print("----------------------------------------")
            if i >= top_candidates:
                break
            else:
                i += 1
    
    def export_homologous_sequences(self, kmer: str, seq1_start: int, seq2_start: int) -> tuple:
        """
        Export the homologous sequences associated with a specific k-mer and start positions.

        Parameters:
        -----------
        kmer : str
            The k-mer for which to export homologous sequences.
        seq1_start : int
            The starting position of the alignment in seq1.
        seq2_start : int
            The starting position of the alignment in seq2.

        Returns:
        --------
        tuple
            A tuple containing the homologous sequences from seq1 and seq2.
        """
        
        for alignment_info in self.alignment_info_list:
            if alignment_info["kmer"] == kmer:
                if alignment_info["seq1_start"] == seq1_start:
                    if alignment_info["seq2_start"] == seq2_start:
                        seq1_end = alignment_info["seq1_end"]
                        seq2_end = alignment_info["seq2_end"]
                        return self.seq1[seq1_start:seq1_end], self.seq2[seq2_start:seq2_end]
        print("Could not locate homologous sequences.")
    
    
    
