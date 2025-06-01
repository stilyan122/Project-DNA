"""
Component 1 – Introduction to Genetic Biology Terms
===================================================

A grab-bag of small, self-contained utilities used in the First Component.
"""

import numpy as np
import random
from collections import Counter

from plotting import quick_line_plot 

def graph_growth_in_cell_division(
    number_of_generations,
    plot_label = "Exponential Growth in Cell Division",
    initial_cells = 1,
    line_color = "y"):
    
    """
    Plot the exponential increase in cell count during mitotic division.
    """

    # Create generations points for the graph
    generations = np.arange(number_of_generations + 1)
    
    # Calculate total cells
    cells = initial_cells * (2 ** generations)

    # Plot the results with the helper function
    quick_line_plot(
        generations,
        cells,
        xlabel="Generations",
        ylabel="Number of Cells",
        title=plot_label,
        label="Cell Division",
        color=line_color,
        xlim=(0, number_of_generations),
        ylim=(0, cells.max() * 1.05),
    )


def compute_dna_entropy(dna_sequence):

    """
    Compute Shannon entropy (in bits) of a single‑stranded DNA sequence.
    """

    # Count occurrences of each nucleotide (A, C, G, T)
    counts = Counter(dna_sequence.upper())

    # Total length of the strand
    seq_len = sum(counts.values())

    # Convert raw counts to probabilities p_i = count / length
    probabilities = (
        count / seq_len for count in counts.values()
    )

    # Shannon entropy (skip p_i = 0 to avoid log2(0))
    entropy_bits = -sum(p * np.log2(p) for p in probabilities if p > 0)

    # Return entropy as function output
    return entropy_bits

def random_dna_sequence(length):
    
    """
    Generate a random DNA sequence of the given length.
    """

    # Build the sequence by randomly choosing one nucleotide per position
    return ''.join(random.choice('ATCG') for _ in range(length))

def plot_entropy_vs_length(start_len, end_len, step):
    
    """
    Plot Shannon entropy of random DNA sequences over a length range.
    """

    # Build the array of sequence lengths
    lengths = np.arange(start_len, end_len + 1, step)

    # Store the entropy value for each length
    entropies = []

    # For every length - generate a random sequence - compute entropy
    for n in lengths:
        seq = random_dna_sequence(n)          
        entropies.append(compute_dna_entropy(seq))

    # Visualise the results with the shared plotting helper
    quick_line_plot(
        lengths,
        entropies,
        xlabel="Sequence Length",
        ylabel="Shannon Entropy",
        title="Entropy vs. Sequence Length (random DNA)",
        label="Random DNA",
        color="b",
        xlim=(0, lengths.max() + 10),
        ylim=(0, max(entropies) + 1),
    )

def mutate_dna_sequence(seq, mutation_rate):
    
    """
    Return a new DNA string after applying random point mutations at the given rate.
    """
    
    # Convert the input string to a NumPy array for easy masking
    dna = np.array(list(seq))

    # Transition <-> transversion lookup tables
    transitions = {
        'A': 'G', 
        'G': 'A', 
        'C': 'T', 
        'T': 'C'
    }
    
    transversions = {
        'A': ['C', 'T'],
        'G': ['C', 'T'],
        'C': ['A', 'G'],
        'T': ['A', 'G'],
    }

    # Decide which positions mutate, then whether each is a transition
    rand_vals = np.random.random(len(dna)) < mutation_rate           # mutate?
    transition_mask = rand_vals & (np.random.random(len(dna)) < 0.5) # half become transitions
    transversion_mask = rand_vals & ~transition_mask                 # the rest become transversions

    # Apply transitions (A <-> G, C <-> T)
    dna[transition_mask] = [transitions[b] for b in dna[transition_mask]]

    # Apply transversions (purine <-> pyrimidine swaps chosen at random)
    dna[transversion_mask] = [
        np.random.choice(transversions[b]) for b in dna[transversion_mask]
    ]

    # Return the mutated sequence as a string
    return ''.join(dna)

def apply_indel_mutations(seq, ins_rate, del_rate):
    
    """
    Create a new DNA string by randomly inserting and deleting bases at the given rates.
    """
    
    # Work with a NumPy array for easy masking
    dna = np.array(list(seq))

    # Deletions: keep bases whose random draw exceeds del_rate
    keep_mask = np.random.random(dna.size) > del_rate
    dna = dna[keep_mask]  # drop deleted positions

    # Insertions: decide where to insert an extra nucleotide (before each kept base)
    insert_mask = np.random.random(dna.size) < ins_rate
    random_bases = np.random.choice(list("ACGT"), insert_mask.sum())  # bases to insert

    # Build the new sequence with insertions interleaved
    output = []
    insert_idx = 0
    for base, do_insert in zip(dna, insert_mask):
        if do_insert:
            output.append(random_bases[insert_idx])  # insert random base first
            insert_idx += 1
        output.append(base)                          # then keep the original

    # Return function output
    return "".join(output)

def simulate_mutation_drift(generations,
                            minimal_DNA_length,
                            point_mutation_rate,
                            insertion_rate,
                            deletion_rate):
    
    """
    Simulate sequence-length drift across generations under point and INDEL mutations.
    """
    
    # Start with a random DNA string of the requested minimal length
    dna = ''.join(np.random.choice(list("ACGT"), size=minimal_DNA_length))
    lengths = [len(dna)]

    # Apply mutations for the desired number of generations
    for _ in range(generations):
        dna = mutate_dna_sequence(dna, point_mutation_rate)          # point mutations
        dna = apply_indel_mutations(dna, insertion_rate, deletion_rate)  # insertions & deletions
        lengths.append(len(dna))                                     # track length drift

    # Visualise length vs. generation with the shared plotting helper
    quick_line_plot(
        range(generations + 1),
        lengths,
        xlabel="Generation",
        ylabel="Sequence Length",
        title="Length drift from point mutations + indels",
        label="DNA length",
        color="b",
    )