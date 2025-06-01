"""
Component 4 – Markov Chains & DNA Mutation Modelling
====================================================

This module will contain the utilities, simulation routines, and analysis
helpers for modelling DNA mutations as Markov chains and other stochastic
processes for the Fourth Component.
"""

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

from plotting import quick_bar_plot 

def jukes_cantor_step(dist, mu):
    
    """
    Apply one Jukes-Cantor Markov-chain step to a nucleotide distribution.
    """
    
    # Build the 4x4 transition matrix
    P = np.full((4, 4), mu / 3)  # off-diagonal entries
    np.fill_diagonal(P, 1 - mu)  # probability of staying the same

    # Treat the input distribution as a row vector and multiply
    v0 = np.asarray(dist, dtype=float)          # [p_A, p_C, p_G, p_T]
    return v0 @ P    

def generate_jc_transition_matrix(mu):
    
    """
    Create a 4x4 Jukes–Cantor transition matrix for a given mutation rate.
    """
    
    # Start with an array filled with off-diagonal probability mu/3
    P = np.full((4, 4), mu / 3)

    # Set diagonal entries to the probability of no mutation (1 − mu)
    np.fill_diagonal(P, 1 - mu)

    return P  # rows sum to 1, order is [A, C, G, T]

def is_stochastic_matrix(P):
    
    """
    Check if each row of matrix P sums to 1 (stochastic matrix).
    """
    
    # Convert input to a NumPy array for reliable summation
    arr = np.asarray(P, dtype=float)

    # Compute the sum of each row
    row_sums = arr.sum(axis=1)

    # Verify every row sum is (approximately) 1
    return bool(np.allclose(row_sums, 1.0, atol=1e-8))

def plot_markov_chain_graph(P, highlight_edges=None, title=""):
    
    """
    Visualise a DNA-mutation Markov chain transition graph.
    """
    
    # Default for highlight_edges if none provided
    if highlight_edges is None:
        highlight_edges = []

    # Define state labels and create a directed multigraph
    state_names = ["A", "C", "G", "T"]
    G = nx.MultiDiGraph()
    G.add_nodes_from(state_names)

    # Add directed edges for each non-self transition with weights
    for i in range(4):
        for j in range(4):
            if i != j:
                G.add_edge(
                    state_names[i],
                    state_names[j],
                    weight=round(P[i][j], 2),
                )

    # Choose a circular layout for clarity
    pos = nx.circular_layout(G)
    plt.figure(figsize=(8, 6))

    # Draw nodes and labels
    nx.draw_networkx_nodes(G, pos, node_color="burlywood", node_size=1000)
    nx.draw_networkx_labels(G, pos, font_size=16)

    # Draw each edge once, with optional highlighting
    drawn = set()
    for u, v, key in G.edges(keys=True):
        if (u, v) in drawn:
            continue
        drawn.add((u, v))

        # Determine if this edge should be highlighted
        is_highlighted = (u, v) in highlight_edges
        edge_color = "red" if is_highlighted else "gray"
        edge_width = 4.0 if is_highlighted else 1.5

        # Compute positions for arrow and label
        dx, dy = pos[v][0] - pos[u][0], pos[v][1] - pos[u][1]
        norm = np.hypot(dx, dy) or 1.0
        dx, dy = dx / norm, dy / norm
        perp_dx, perp_dy = -dy, dx

        # Midpoint offset so parallel edges won’t overlap (not needed for a single arrow)
        control_x = (pos[u][0] + pos[v][0]) / 2
        control_y = (pos[u][1] + pos[v][1]) / 2

        # Draw the directed edge with an arrow
        arrowprops = dict(
            arrowstyle="-|>",
            color=edge_color,
            lw=edge_width,
            shrinkA=15,
            shrinkB=15,
        )
        plt.annotate(
            "",
            xy=pos[v],
            xytext=pos[u],
            arrowprops=arrowprops,
            annotation_clip=False,
            xycoords="data",
            textcoords="data",
        )

        # Place the weight label near the midpoint
        weight_label = str(G[u][v][key]["weight"])
        plt.text(
            control_x,
            control_y,
            weight_label,
            fontsize=10,
            ha="center",
            va="center",
            backgroundcolor="white",
            color="black",
        )

    # Show self-transition probabilities above each node
    for i, name in enumerate(state_names):
        x, y = pos[name]
        plt.text(x, y + 0.10, round(P[i][i], 2), fontsize=11, ha="center", color="black")

    # Add title and hide axes
    plt.text(
        0.5,
        1.05,
        title,
        fontsize=14,
        ha="center",
        va="bottom",
        transform=plt.gca().transAxes,
    )
    
    plt.axis("off")
    plt.tight_layout()
    plt.show()

def highlight_transition(P, from_idx, to_idx):
    
    """
    Highlight a single transition in the mutation Markov chain graph.
    """
    
    # Map index to nucleotide label
    state_names = ["A", "C", "G", "T"]
    src = state_names[from_idx]
    tgt = state_names[to_idx]
    edge = (src, tgt)

    # Reuse the Markov-chain visualiser with this one edge emphasized
    plot_markov_chain_graph(
        P,
        highlight_edges=[edge],
        title=rf"Accessibility: {src} $\rightarrow$ {tgt}",
    )

def highlight_communication(P, from_idx, to_idx):
    """
    Highlight bidirectional communication between two states in the Markov chain graph.
    """
    
    # Map numeric indices to nucleotide labels
    state_names = ["A", "C", "G", "T"]
    src = state_names[from_idx]
    tgt = state_names[to_idx]

    # Define the two directed edges to highlight (src→tgt and tgt→src)
    edge_forward = (src, tgt)
    edge_backward = (tgt, src)

    # Reuse the graph-plotting helper to emphasize both edges
    plot_markov_chain_graph(
        P,
        highlight_edges=[edge_forward, edge_backward],
        title=rf"Communication: {src} $\leftrightarrow$ {tgt}",
    )

def highlight_irreducibility(P):
    
    """
    Highlight all possible transitions to demonstrate irreducibility of the Markov chain.
    """
    
    # Define state labels and list all directed edges between distinct states
    state_names = ["A", "C", "G", "T"]
    all_edges = [(a, b) for a in state_names for b in state_names if a != b]

    # Reuse the graph-plotting helper to emphasize every edge
    plot_markov_chain_graph(
        P,
        highlight_edges=all_edges,
        title="Irreducibility"
    )

def find_stationary_distribution(P):
    
    """
    Compute the stationary distribution of a transition matrix P.
    """
    
    arr = np.array(P, dtype=float)

    # If P is identity, there's no unique non-trivial stationary distribution
    if np.allclose(arr, np.eye(arr.shape[0])):
        return None

    # Set up the linear system: pi * P = pi and sum(pi) = 1
    n = arr.shape[0]
    A = arr.T - np.eye(n)        # (P^T − I)
    A[-1] = np.ones(n)           # replace last row by [1, 1, ..., 1]
    b = np.zeros(n)
    b[-1] = 1                    # enforces sum(pi) = 1

    # Solve for pi
    pi = np.linalg.solve(A, b)

    # Return function output
    return pi

def plot_stationary_distribution(pi):
    
    """
    Plot the stationary distribution pi as a bar chart.
    """
    
    # Define state labels
    state_names = ["A", "C", "G", "T"]

    # Draw bars for pi using the shared helper (keep figure open for annotations)
    quick_bar_plot(
        categories=state_names,
        values=pi.tolist(),
        title=r"Stationary Distribution $\pi$",
        ylabel="Probability",
        colors=["cornflowerblue"] * 4,
        ylim=(0, 1),
        grid=True,
        show=False,          # add text labels next
    )

    # Annotate each bar with its numeric probability
    bars = plt.gca().patches
    for bar, prob in zip(bars, pi):
        plt.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.02,
            f"{prob:.2f}",
            ha="center",
            fontsize=10,
        )

    # Finalize and display
    plt.title(r"Stationary Distribution $\pi$", pad=15)
    plt.tight_layout()
    plt.show()

def compute_and_plot_stationary(P):
    
    """
    Compute and display the stationary distribution of a transition matrix.
    """
    
    # Find the stationary distribution (pi) using the solver
    pi = find_stationary_distribution(P)

    # If the matrix is not identity, plot pi; otherwise print a message
    if pi is not None:
        # Use the bar-plot helper for a clean visualization
        state_names = ["A", "C", "G", "T"]
        quick_bar_plot(
            categories=state_names,
            values=pi.tolist(),
            title=r"Stationary Distribution $\pi$",
            ylabel="Probability",
            colors=["cornflowerblue"] * 4,
            ylim=(0, 1),
            grid=True,
            show=False,  # keep open to annotate
        )

        # Annotate each bar with its probability
        ax = plt.gca()
        for bar, prob in zip(ax.patches, pi):
            plt.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.02,
                f"{prob:.2f}",
                ha="center",
                fontsize=10,
            )

        plt.tight_layout()
        plt.show()
    else:
        print("Mutation rate is zero — the matrix is identity. Any initial distribution is stationary.")