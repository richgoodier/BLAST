import numpy as np
import matplotlib.pyplot as plt

def plot_alignment(seq: str, matches: np.ndarray, title: str, step: int = 100) -> None:
    """
    Plot the alignment between two sequences based on match scores.

    Parameters:
    seq (str): The first sequence string.
    matches (np.ndarray): Match scores array from alignment.
    title (str): Title of the plot.
    step (int): Step size for downsampling the data points in the plot.

    Returns:
    None
    """
    plt.figure(figsize=(10, 6))
    alignment_shift = [i - len(seq) for i in range(len(matches))]
    plt.plot(alignment_shift[::step], matches[::step])
    plt.title(title)
    plt.xlabel('Alignment Shift')
    plt.ylabel('Number of Nucleotides in Alignment')
    plt.grid(True)
    plt.show()


def plot_frequency(freq: np.ndarray, title: str) -> None:
    """
    Plot the frequency of matches for a subset of alignment shifts.

    Parameters:
    freq (np.ndarray): Frequency of matches array.
    title (str): Title of the plot.

    Returns:
    None
    """
    clip = int(len(freq) * 0.1)
    max_alignment = np.argmax(freq[clip:-clip]) + clip # Remove wings to avoid edge artifacts

    breadth = int(0.35 * len(freq)) # Zoom in on the max_alignment and a bit on either side.
    step = max(1, int(breadth / 100)) # Reduce the number of points plotted

    start = max(0, max_alignment - breadth) # Define the range of indices around max_alignment
    end = min(len(freq), max_alignment + breadth)

    plt.figure(figsize=(10, 6))
    alignment_shift = [i - len(freq) // 2 for i in range(len(freq))]
    plt.plot(alignment_shift[start:end:step], freq[start:end:step])
    plt.title(title)
    plt.xlabel('Alignment Shift')
    plt.ylabel('Percentage of Nucleotides in Alignment')
    plt.grid(True)
    plt.show()