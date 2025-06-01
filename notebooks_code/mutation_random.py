"""
Component 3 – Mathematical Representations of Biological Structures
===================================================================

This module will gather all helper functions, classes, and utilities used in the Second Component.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import binom

from plotting import quick_line_plot, quick_bar_plot 

def demo_die_probability(event, n_rolls):
    
    """
    Demonstrate theoretical vs. empirical probability when rolling a fair die.
    """
    
    # Define the sample space (faces 1-6) and its theoretical probability
    faces = np.arange(1, 7)           # sample space omega = {1,2,3,4,5,6}
    p_theoretical = 1 / len(faces)    # p(face) = 1/6 for a fair die

    # Simulate n_rolls independent rolls
    rolls = np.random.choice(faces, size=n_rolls, replace=True)
    p_empirical = np.mean(rolls == event)   # empirical probability of the target face

    # Print the comparison
    print(f"Theoretical P(roll = {event}): {p_theoretical:.8f}")
    print(f"Empirical P(roll = {event}): {p_empirical:.8f}")

    # Visualise the counts for each face
    counts = [np.sum(rolls == face) for face in faces]

    quick_bar_plot(
        categories=faces,
        values=counts,
        xlabel="Die face",
        ylabel="Count",
        title=f"Empirical Frequencies of Die Rolls ({n_rolls} trials)",
        colors=["skyblue"] * 6,
        grid=True,
        alpha=0.8,
    )

def simulate_bernoulli_mutations(mu, n_trials):
    
    """
    Simulate Bernoulli mutation trials (1 = mutation, 0 = no mutation) and plot the results.
    """
    
    # Generate n_trials Bernoulli outcomes
    samples = (np.random.rand(n_trials) < mu).astype(int)   # 0 / 1 array

    # Count how many 0s and 1s occurred
    counts = np.bincount(samples, minlength=2)              # [n_no_mut, n_mut]

    # Visualise with the shared bar-plot helper
    quick_bar_plot(
        categories=["No mutation (0)", "Mutation (1)"],
        values=counts.tolist(),
        xlabel="Outcome",
        ylabel="Count",
        title=f"Bernoulli Trials (μ = {mu}, n = {n_trials})",
        colors=["gray", "blue"],
        alpha=0.7,
        grid=True,
    )

def simulate_binomial_mutations(n, mu, n_sims):
    
    """
    Compare the theoretical Binomial (n, mu) PMF with simulated mutation counts.
    """
    
    # Theoretical PMF
    k = np.arange(n + 1)          # possible mutation counts 0 … n
    pmf = binom.pmf(k, n, mu)       # Binomial probability-mass function

    # Simulation: draw n_sims samples of "mutations out of n trials"
    samples = np.random.binomial(n, mu, size=n_sims)

    # Plot: line for PMF, histogram for simulation
    plt.figure(figsize=(8, 4))
    plt.plot(k, pmf, "o-", color="red", label="Theoretical PMF")
    plt.hist(
        samples,
        bins=np.arange(-0.5, n + 1.5),
        density=True,
        alpha=0.5,
        color="blue",
        label="Simulated",
    )

    plt.xlabel("Total mutations (k)")
    plt.ylabel("Probability")
    plt.title(f"Binomial Distribution (n = {n}, μ = {mu})")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

def simulate_mutation_type_distribution(p_transition, n_bases):
    
    """
    Simulate transition vs transversion outcomes and show their empirical frequencies.
    """
    
    # Draw n_bases uniform random numbers; classify each as transition (1) or transversion (0)
    is_transition = (np.random.rand(n_bases) < p_transition).astype(int)

    # Compute empirical frequencies
    counts = np.bincount(is_transition, minlength=2)     # [transversion, transition]
    freqs  = counts / n_bases

    # Visualise with the shared bar-plot helper
    quick_bar_plot(
        categories=["Transversion", "Transition"],
        values=freqs.tolist(),
        ylabel="Empirical Frequency",
        title=f"Categorical Distribution (p_transition = {p_transition}, n = {n_bases})",
        colors=["orange", "green"],
        ylim=(0, 1.05),
        alpha=0.8,
        grid=True,
        show=False
    )

    # Annotate bars with precise frequencies
    for x, f in enumerate(freqs):
        plt.text(x, f + 0.01, f"{f:.5f}", ha="center")
    plt.tight_layout()
    plt.show()