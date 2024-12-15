from Bio import SeqIO

def fasta_to_str(filepath):
    for record in SeqIO.parse(filepath, 'fasta'):
        return record.seq

def rev_complement(seq):
    complement = {"G": "C",
                  "C": "G",
                  "A": "T",
                  "T": "A",
                  "N": "N",
                  "g": "c",
                  "c": "g",
                  "a": "t",
                  "t": "a"}
    seq_len = len(seq)
    output = ""
    
    for i in range(seq_len):
        output += complement[seq[seq_len-i-1]]
    
    return output

def search_for_kmer(seq1, seq2, kmer):
    def find_kmer_indices(sequence, kmer):
        """Find all indices of kmer in a given sequence."""
        kmer_length = len(kmer)
        return [i for i in range(len(sequence) - kmer_length + 1) if sequence[i:i + kmer_length] == kmer]

    # Find indices in both sequences
    seq1_indices = find_kmer_indices(seq1, kmer)
    seq2_indices = find_kmer_indices(seq2, kmer)

    # Store results in dictionaries
    seq1_index = {kmer: seq1_indices}
    seq2_index = {kmer: seq2_indices}

    # Print summary
    print(f"{kmer} ---- {len(seq1_indices)} ---- {len(seq2_indices)}")

    return seq1_index, seq2_index


def check_homologous_sequences(seq1, seq2):
    if len(seq1) != len(seq2):
        print("Sequences must be of equal lengths")
        return
    
    length_seq = len(seq1)
    num_matches = 0
    for i in range(length_seq):
        if seq1[i] == seq2[i]:
            num_matches += 1
    
    print(f"A total of {num_matches} of {len(seq1)} bases in alignment.")
    print(f"Fraction aligned: {num_matches / length_seq}")
