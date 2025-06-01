"""
Component 2 – Mutation as a Random Process
==========================================

This module will collect the utilities and simulators used in the Third Component
of the project, where we treat mutation as a stochastic process (Markov
chains).
"""

import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

from plotting import quick_line_plot, quick_bar_plot 

from bio_intro import mutate_dna_sequence

def encode_sequence(seq, mapping):
    
    """
    Convert a DNA/RNA string into a NumPy array of mapped values.
    """
    
    # Build a list using the provided mapping, then cast to NumPy for vector ops
    return np.array([mapping[b] for b in seq])

def visualize_sequence_space(max_len):
    
    """
    Plot the explosive growth of possible DNA/RNA sequences as length increases.
    """
    
    #  Build length vector (1 - max_len) and compute 4^n possibilities
    n = np.arange(1, max_len + 1)
    counts = 4.0 ** n

    # Draw the line with the shared helper (keep figure open)
    quick_line_plot(
        n,
        counts,
        xlabel="Sequence length n",
        ylabel="Number of possible sequences",
        title="DNA / RNA Sequence-Space Explosion",
        label="4ⁿ",
        color="b",
        show=False,           # we still need to tweak the axes
    )

    # Log-scale y-axis and finalise the plot
    plt.yscale("log")
    plt.legend()
    plt.show()

def encode_codon_triplets(sequence, mapping):
    
    """
    Transform a DNA/RNA string into a list of 3-tuple codon encodings.
    """
    
    # First map each nucleotide to its numeric code
    encoded = encode_sequence(sequence, mapping)     # we rely on helper defined earlier

    codons = []  # collect (b1, b2, b3) triplets

    # Walk the encoded array three bases at a time, skip leftovers
    for i in range(0, len(encoded) - len(encoded) % 3, 3):
        triple = encoded[i : i + 3]
        codons.append((int(triple[0]), int(triple[1]), int(triple[2])))

    # Return function output
    return codons

def encode_protein_sequence(seq, mapping):
    
    """
    Convert an amino-acid string into a NumPy array of mapped values.
    """
    
    # Map each residue's one-letter code to its numeric representation
    return np.array([mapping[a] for a in seq])

def translate_protein_to_dna(protein_seq, codon_map):
    
    """
    Reverse-translate a protein string into a DNA sequence using a codon lookup.
    """
    
    # Build a list of codons corresponding to each amino acid
    dna_codons = []
    for aa in protein_seq:
        codon = codon_map.get(aa, "NNN")   # fallback to 'NNN' if AA not in map
        dna_codons.append(codon)

    # Join codons with single spaces for readability (e.g. "ATG GCT ...")
    return " ".join(dna_codons)

def nucleotide_frequencies(seq):
    
    """
    Return the relative frequency of each nucleotide in a DNA string.
    """
    
    # Count how many times each base appears (case-insensitive)
    counts = Counter(seq.upper())

    # Convert counts to fractions of the total length
    total = len(seq)
    return {b: counts.get(b, 0) / total for b in "ACGT"}

def simulate_mutation_frequency_history(initial_seq, mutation_rate, generations):
   
    """
    Record nucleotide-frequency dynamics across generations under point mutation.
    """
    
    # Initialise the record with the base composition of the starting strand
    history = []
    seq = initial_seq

    # For each generation - store frequencies - mutate the sequence
    for _ in range(generations + 1):
        history.append(nucleotide_frequencies(seq))          # current composition
        seq = mutate_dna_sequence(seq, mutation_rate)        # apply point mutations

    # Return the list of frequency dictionaries (one entry per generation)
    return history

def plot_mutation_dynamics(history):
    
    """
    Visualise nucleotide-frequency trajectories and the final composition histogram.
    """
    
    # Line-plot: frequency of each base across generations
    generations = list(range(len(history)))
    for base in "ACGT":
        freqs = [h[base] for h in history]          # pull base column
        quick_line_plot(
            generations,
            freqs,
            xlabel="Generation",
            ylabel="Frequency",
            title="Point-mutation dynamics by base",
            label=base,
            grid=True,
            show=False,                             # overlay all four lines first
        )
        
    plt.legend()
    plt.show()

    # Bar-chart: final generation's nucleotide composition
    
    final_freqs = history[-1]
    quick_bar_plot(
        categories=list("ACGT"),
        values=[final_freqs[b] for b in "ACGT"],
        xlabel="Base",
        ylabel="Final Frequency",
        title="Final nucleotide frequencies",
        colors=["b", "g", "r", "y"],
        grid=True,
    )

def mutate_sequence_matrix(seq, matrix):
    
    """
    Mutate a DNA string according to a 4 x 4 substitution-probability matrix.
    """
   
    mapping = {"A": 0, "C": 1, "G": 2, "T": 3}   # base -> numeric index
    bases   = ["A", "C", "G", "T"]               # index -> base

    # Encode the input sequence as numeric indices
    encoded = encode_sequence(seq, mapping)       

    # Sample a new base for every position using the row in the matrix
    mutated = np.array([np.random.choice(4, p=matrix[i]) for i in encoded])

    # Decode back to characters and return the mutated DNA
    return "".join(bases[i] for i in mutated)

def simulate_matrix_mutation_history(seq, matrix, generations):
    
    """
    Track nucleotide-frequency history for matrix-driven point mutations.
    """
    
    history = []             # will hold one dict per generation
    dna = seq                # work on a local copy

    for _ in range(generations + 1):
        history.append(nucleotide_frequencies(dna))         # record composition
        dna = mutate_sequence_matrix(dna, matrix)           # apply one round of mutations

    # Return the desired output
    return history  

def compare_mutation_histories(history_flip, history_matrix):
    
    """
    Compare base-frequency trajectories for flip- and matrix-mutation models.
    """

    bases  = "ACGT"
    colors = ["b", "g", "r", "y"]
    gens   = np.arange(len(history_flip))

    # Trajectory plot
    plt.figure(figsize=(10, 5))
    for base, col in zip(bases, colors):
        # Flip model - solid line
        plt.plot(
            gens,
            [h[base] for h in history_flip],
            marker="o",
            linestyle="-",
            color=col,
            label=f"{base} (flip)",
        )
        
        # Matrix model - dashed line
        plt.plot(
            gens,
            [h[base] for h in history_matrix],
            marker="s",
            linestyle="--",
            color=col,
            label=f"{base} (matrix)",
        )

    plt.xlabel("Generation")
    plt.ylabel("Frequency")
    plt.title("Flip vs. Matrix Mutation Dynamics")
    plt.grid(True)
    plt.legend(ncol=2, fontsize="small")
    plt.tight_layout()
    plt.show()

    # Final-generation histogram 
    final_flip = history_flip[-1]
    final_mat  = history_matrix[-1]

    bases_list = list(bases)
    x          = np.arange(len(bases_list))
    width      = 0.35                       # equal width for both series

    plt.figure(figsize=(6, 4))
    plt.bar(
        x - width / 2,
        [final_flip[b] for b in bases_list],
        width=width,
        color="C0",
        alpha=0.7,
        label="Flip",
    )
    plt.bar(
        x + width / 2,
        [final_mat[b] for b in bases_list],
        width=width,
        color="C1",
        alpha=0.9,
        label="Matrix",
    )

    plt.xticks(x, bases_list)
    yticks = np.arange(0, 1.01, 0.05)
    plt.yticks(yticks, [f"{t:.2f}" for t in yticks])
    plt.ylim(0, 1)
    plt.xlabel("Base")
    plt.ylabel("Final Frequency")
    plt.title("Final Nucleotide Frequencies: Flip vs. Matrix")
    plt.grid(axis="y")
    plt.legend()
    plt.tight_layout()
    plt.show()