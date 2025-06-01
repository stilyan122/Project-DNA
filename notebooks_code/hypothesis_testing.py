"""
Component 5 – Hypothesis Testing
================================

This module will house the utilities, statistical tests, and analysis functions
for biological hypothesis testing stated in the Fifth Component
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chisquare

from plotting import quick_bar_plot, quick_line_plot
from markov_mutations import generate_jc_transition_matrix, find_stationary_distribution
from bio_intro import random_dna_sequence
from bio_structures import simulate_matrix_mutation_history

def simulate_mutations_and_stationary(mu, seq_length, generations, seed=None):
    
    """
    Simulate JC69 mutations to compute the stationary distribution and frequency history.
    """
    
    # Optionally set random seed for reproducibility
    if seed is not None:
        np.random.seed(seed)

    # Build the Jukes–Cantor transition matrix for mutation rate mu
    P = generate_jc_transition_matrix(mu)

    # Compute the stationary distribution pi of P
    pi = find_stationary_distribution(P)
    
    if pi is None:
        # If P is identity (mu = 0), the uniform distribution is stationary
        pi = np.array([0.25, 0.25, 0.25, 0.25])

    # Generate a random DNA string of length N and simulate mutations over generations
    initial_seq = random_dna_sequence(seq_length)
    freqs_history = simulate_matrix_mutation_history(initial_seq, P, generations)

    # Return the stationary distribution and the full frequency history
    return pi, freqs_history

def compute_l1_distance_history(freqs_history, pi):
    
    """
    Compute the L1 distance from the stationary distribution at each generation.
    """
    
    # Define bases and map pi to a dictionary for quick lookup
    bases   = ["A", "C", "G", "T"]
    pi_dict = dict(zip(bases, pi))

    # For each generation's frequency dict, sum absolute differences |f(b) − pi(b)|
    return [
        sum(abs(freqs[b] - pi_dict[b]) for b in bases)
        for freqs in freqs_history
    ]

def plot_l1_distances(l1_history, epsilon):
    
    """
    Plot the L1-distance curve versus generation with an epsilon threshold line.
    """
    
    # Draw the L1-distance curve using the shared line-plot helper
    quick_line_plot(
        range(len(l1_history)),
        l1_history,
        xlabel="Generation $t$",
        ylabel="L₁ Distance",
        title="L₁ Distances to Stationary Distribution",
        label=r'$\|\hat\pi_t - \pi\|_1$',
        color="b",
        grid=True,
        show=False,               # overlay the epislon line next
    )

    # Draw a horizontal epsilon line at the specified threshold
    plt.axhline(y=epsilon, linestyle="--", color="r", label=f"ε = {epsilon}")

    # Finalize labels, legend, and display
    plt.legend()
    plt.tight_layout()
    plt.show()

def print_frequency_comparison(final_freqs, pi, bases):
    
    """
    Print a table comparing observed vs. theoretical nucleotide frequencies.
    """
    
    # Print header row
    print("Base\tObserved\tTheoretical\t|Δ|")
    # For each base, show observed frequency, pi, and absolute difference
    for b in bases:
        obs  = final_freqs[b]
        theo = pi[bases.index(b)]
        diff = abs(obs - theo)
        print(f"{b}\t{obs:.4f}\t\t{theo:.2f}\t\t{diff:.4f}")

def plot_frequency_comparison(final_freqs, pi, bases):
    
    """
    Draw a grouped bar chart comparing observed vs. theoretical base frequencies.
    """
    
    # Extract observed frequencies in the same order as bases
    observed_vals = [final_freqs[b] for b in bases]
    theoretical_vals = list(pi)

    # Define positions and bar width for two groups
    x = np.arange(len(bases))
    width = 0.35

    # Create the figure and plot both series side by side
    plt.figure(figsize=(6, 4))
    
    # Theoretical frequencies on the left side of each group
    plt.bar(
        x - width / 2,
        theoretical_vals,
        width=width,
        label="Theoretical",
        alpha=0.6,
        color="cornflowerblue",
        edgecolor="black",
    )
    # Observed frequencies on the right side of each group
    plt.bar(
        x + width / 2,
        observed_vals,
        width=width,
        label="Observed",
        alpha=0.6,
        color="orange",
        edgecolor="black",
    )

    # Finalize axis labels, ticks, legend, and layout
    plt.xticks(x, bases)
    plt.ylabel("Frequency")
    plt.title("Observed vs. Theoretical Base Frequencies")
    plt.legend()
    plt.ylim(0, 1)
    plt.grid(axis="y", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.show()

def run_mutation_analysis(mu, seq_length, generations, epsilon, seed, bases):
    
    """
    Execute the JC69 mutation simulation, compute L1-distance history, and display summaries.
    """
    
    # Simulate mutations and get stationary distribution + frequency history
    pi, freqs_history = simulate_mutations_and_stationary(mu, seq_length, generations, seed)

    # Compute the L1-distance at each generation from the stationary distribution
    l1_history = compute_l1_distance_history(freqs_history, pi)

    # Plot the L1-distance curve with the epsilon threshold line
    plot_l1_distances(l1_history, epsilon)

    # Extract final-generation frequencies
    final_freqs = freqs_history[-1]

    # Print a tabular summary comparing observed vs. theoretical frequencies
    print_frequency_comparison(final_freqs, pi, bases)

    # Draw a grouped bar chart of observed vs. theoretical frequencies
    plot_frequency_comparison(final_freqs, pi, bases)

def compare_multiple_mutation_rates(seq_length, generations, epsilon):
    
    """
    Run the analysis pipeline for several preset mutation rates and display results.
    """
    
    # Define the nucleotide order
    bases = ["A", "C", "G", "T"]

    # Loop over a set of mutation rates and run the full analysis for each
    for mu in [0.0, 0.25, 0.5, 0.75, 1.0]:
        print(f"\n=== Mutation rate mu = {mu} ===")
        run_mutation_analysis(mu, seq_length, generations, epsilon, seed=42, bases=bases)

def perform_chi_square_test(final_freqs, seq_length):
    
    """
    Perform a chi-square goodness-of-fit test comparing observed vs. uniform frequencies.
    """
    
    # Build the observed counts array in A, C, G, T order
    observed = np.array([final_freqs[b] * seq_length for b in ["A", "C", "G", "T"]])

    # The expected counts under uniform (stationary) distribution
    expected = np.full(4, seq_length / 4)

    # Run SciPy's chi-square test
    chi2_stat, p_value = chisquare(observed, f_exp=expected)

    # Return desired output
    return chi2_stat, p_value

def run_chi_square_experiments(mu_list, seq_length, generations, alpha=0.05, seed=None):
    
    """
    For each mu in mu_list, simulate JC69 mutations, run χ² on final counts, and report results.
    """
    
    # Optionally set random seed for reproducibility
    if seed is not None:
        np.random.seed(seed)

    # Define bases for indexing-and total N for count conversion
    bases = ["A", "C", "G", "T"]
    N = seq_length

    # Loop over each mutation rate mu
    for mu in mu_list:
        print(f"\n=== Mutation rate mu = {mu:.2f} ===")

        # Build JC transition matrix and find stationary distribution π
        P  = generate_jc_transition_matrix(mu)
        pi = find_stationary_distribution(P)
        if pi is None:
            # If mu = 0 (identity matrix), uniform π = [1/4, 1/4, 1/4, 1/4]
            pi = np.array([0.25] * 4)

        # Generate a random DNA sequence and simulate mutations
        initial_seq     = random_dna_sequence(N)
        freqs_history = simulate_matrix_mutation_history(initial_seq, P, generations)

        # Extract final-generation frequencies and perform χ² test
        final_freqs = freqs_history[-1]
        chi2_stat, p_value = perform_chi_square_test(final_freqs, N)

        # Print results and decision
        
        print(f"χ² = {chi2_stat:.2f}, p = {p_value:.4f}")
        
        if p_value < alpha:
            print("Reject H0: final counts differ significantly from uniform.")
        else:
            print("Do not reject H0: final counts are consistent with uniform.")

def compute_kl_divergence(freqs, pi):
    
    """
    Compute KL divergence in bits.
    """
    
    bases = ["A", "C", "G", "T"]
    kl = 0.0

    # Sum over each base: p_obs * log2(p_obs / p_th), skipping zeroes
    for i, b in enumerate(bases):
        p_obs = freqs[b]
        p_th  = pi[i]
        if p_obs > 0:
            kl += p_obs * np.log2(p_obs / p_th)

    # Return function output
    return kl

def simulate_and_plot_kl_divergence(mu, seq_length, generations, seed=None):
    
    """
    Simulate JC69 mutations, compute KL divergence to pi each generation, and plot the curve.
    """
    
    # Optionally set random seed for reproducibility
    if seed is not None:
        np.random.seed(seed)

    # Build Jukes–Cantor transition matrix and find its stationary distribution π
    P = generate_jc_transition_matrix(mu)
    pi = find_stationary_distribution(P)
    if pi is None:
        # If mu = 0 (identity matrix), uniform π = [1/4, 1/4, 1/4, 1/4]
        pi = np.array([0.25, 0.25, 0.25, 0.25])

    # Generate initial random DNA sequence and simulate mutation-frequency history
    initial_seq    = random_dna_sequence(seq_length)
    freqs_history = simulate_matrix_mutation_history(initial_seq, P, generations)

    # Compute KL divergence at each generation
    kl_history = [compute_kl_divergence(freqs, pi) for freqs in freqs_history]

    # Plot KL divergence vs. generation using shared line-plot helper
    quick_line_plot(
        range(len(kl_history)),
        kl_history,
        xlabel="Generation $t$",
        ylabel="KL Divergence (bits)",
        title="KL Divergence to Uniform Stationary Distribution",
        label=r"$D_{\mathrm{KL}}(\hat{\pi}_t \,\|\, \pi)$",
        color="purple",
        grid=True,
    )

def bootstrap_frequency_ci(final_states, n_boot=1000, ci=95):
    
    """
    Compute bootstrap confidence intervals for nucleotide frequencies.
    """
    
    # Number of sequences and base labels
    N = len(final_states)
    bases = ["A", "C", "G", "T"]

    # Preallocate array to store bootstrap frequency estimates (shape: n_boot x 4)
    boot_freqs = np.zeros((n_boot, 4), dtype=float)

    # Perform n_boot bootstrap replicates
    for i in range(n_boot):
        # Resample with replacement from final_states (integers 0–3)
        sample = np.random.choice(final_states, size=N, replace=True)
        
        # Count occurrences of each base code 
        counts = np.bincount(sample, minlength=4)
        
        # Convert counts to frequencies and store
        boot_freqs[i] = counts / N

    # Determine percentile cutoffs for the desired confidence interval
    lower_pct = (100 - ci) / 2
    upper_pct = 100 - lower_pct

    # Compute CI bounds for each base
    ci_dict = {}
    for j, b in enumerate(bases):
        lower = np.percentile(boot_freqs[:, j], lower_pct)
        upper = np.percentile(boot_freqs[:, j], upper_pct)
        ci_dict[b] = (lower, upper)

    # Return desired output
    return ci_dict

def run_bootstrap_ci(mu, seq_length, generations, n_boot=1000, ci=95, seed=None):
    
    """
    Simulate JC69 mutations, then bootstrap CIs for final base frequencies.
    """
    
    # Optionally set the random seed
    if seed is not None:
        np.random.seed(seed)

    # Build Jukes–Cantor transition matrix and find its stationary distribution pi
    P = generate_jc_transition_matrix(mu)
    pi = find_stationary_distribution(P)
    if pi is None:
        # If mu = 0 (identity matrix), the stationary distribution is uniform
        pi = np.array([0.25, 0.25, 0.25, 0.25])

    # Generate an initial random DNA sequence as integer codes 0–3
    initial_states = np.random.choice(4, size=seq_length)

    # Simulate the mutation process over the given number of generations
    history = []  # will hold one array of codes per generation
    current = initial_states.copy()
    for _ in range(generations + 1):
        history.append(current.copy())
        # For each position, sample a new base code according to the row in P
        current = np.array([np.random.choice(4, p=P[code]) for code in current])

    # Extract final-generation array of integer codes
    final_states = history[-1]  # shape = (seq_length,)

    # Compute bootstrap CIs for each base frequency
    ci_dict = bootstrap_frequency_ci(final_states, n_boot=n_boot, ci=ci)

    # Print the CIs alongside the theoretical π
    print("Bootstrap Confidence Intervals for Final Base Frequencies:")
    for j, base in enumerate(["A", "C", "G", "T"]):
        low, high = ci_dict[base]
        print(f"  {base}: [{low:.3f}, {high:.3f}] (π = {pi[j]:.3f})")

    # Show point estimates of final frequencies
    counts = np.bincount(final_states, minlength=4)
    final_freqs = counts / seq_length
    print("\nPoint Estimates of Final Frequencies:")
    for j, base in enumerate(["A", "C", "G", "T"]):
        print(f"  {base}: {final_freqs[j]:.4f}")

def get_final_states_from_freqs(final_freqs, N):
    
    """
    Convert a final-frequency dictionary into an integer-coded array of length N.
    """
    
    bases = ["A", "C", "G", "T"]

    # Compute raw integer counts by rounding freq * N
    raw_counts = {b: int(round(final_freqs[b] * N)) for b in bases}
    total = sum(raw_counts.values())

    # Adjust counts if rounding error causes a mismatch with N
    if total != N:
        # Compute fractional parts for each base’s ideal count
        frac_parts = {
            b: (final_freqs[b] * N) - np.floor(final_freqs[b] * N)
            for b in bases
        }
        
        # Sort bases by descending fractional part
        sorted_bases = sorted(bases, key=lambda b: frac_parts[b], reverse=True)
        diff = N - total
        idx = 0
        
        # Distribute the rounding error (+1 or -1) based on largest fractional parts
        while total != N:
            b = sorted_bases[idx % 4]
            if diff > 0:
                raw_counts[b] += 1
                total += 1
            else:
                if raw_counts[b] > 0:
                    raw_counts[b] -= 1
                    total -= 1
            idx += 1

    # Build the final_states array: repeat each code according to count
    final_states = np.concatenate([
        np.full(raw_counts[b], code, dtype=np.int8)
        for code, b in enumerate(bases)
    ])

    # Shuffle so the integer-coded sequence is randomized
    np.random.shuffle(final_states)
    return final_states

def run_bootstrap_ci_pipeline(mu, seq_length, generations, n_boot=1000, ci=95, seed=None):
    
    """
    Simulate JC69 mutations, convert final frequencies to int‐codes, and bootstrap CI.
    """
    
    # Optionally set random seed for reproducibility
    if seed is not None:
        np.random.seed(seed)

    # Build the Jukes–Cantor transition matrix and compute stationary pi
    P  = generate_jc_transition_matrix(mu)
    pi = find_stationary_distribution(P)
    if pi is None:
        pi = np.array([0.25, 0.25, 0.25, 0.25])

    # Generate an initial random DNA sequence and simulate mutation frequencies
    initial_seq = random_dna_sequence(seq_length)
    freqs_history = simulate_matrix_mutation_history(initial_seq, P, generations)

    # Extract the final‐generation frequency dictionary
    final_freqs = freqs_history[-1]  # {'A':…, 'C':…, 'G':…, 'T':…}

    # Convert final_freqs to an array of integer codes
    final_states = get_final_states_from_freqs(final_freqs, seq_length)

    # Compute bootstrap confidence intervals for each base frequency
    ci_dict = bootstrap_frequency_ci(final_states, n_boot=n_boot, ci=ci)

    # Print the CI results alongside the theoretical π
    print(f"Bootstrap {ci}% Confidence Intervals at Generation {generations}:")
    for j, base in enumerate(["A", "C", "G", "T"]):
        low, high = ci_dict[base]
        print(f"  {base}: [{low:.3f}, {high:.3f}] (π = {pi[j]:.3f})")

    # Also display point estimates of the final frequencies
    counts = np.bincount(final_states, minlength=4)
    estimated = counts / seq_length
    print("\nPoint Estimates of Final Frequencies:")
    for j, base in enumerate(["A", "C", "G", "T"]):
        print(f"  {base}: {estimated[j]:.4f}")