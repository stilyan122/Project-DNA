"""
Helper plotting functions for the whole project.
"""

import matplotlib.pyplot as plt
from typing import Iterable, Sequence, Tuple, List

def quick_line_plot(
    x: Sequence[float] | Iterable[float],
    y: Sequence[float] | Iterable[float],
    *,
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
    label: str | None = None,
    color: str | None = None,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    show: bool = True,
    grid: bool = True
):
    
    """
    Common line‑plot helper used throughout the project.
    """
    
    plt.plot(x, y, label=label, color=color)
    if xlabel:
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)
    if xlim:
        plt.xlim(*xlim)
    if ylim:
        plt.ylim(*ylim)
    if title:
        plt.title(title)
    if label:
        plt.legend()
    if grid:
        plt.grid(True)
    plt.tight_layout()
    if show:
        plt.show()

def quick_bar_plot(
    categories: Sequence[str],
    values: Sequence[float],
    *,
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
    colors: List[str] | None = None,
    ylim: Tuple[float, float] | None = None,
    alpha: float = 0.8,
    show: bool = True,
    grid: bool = True,
):
    
    """
    Minimal wrapper for a single‑series bar chart.
    """

    plt.figure(figsize=(6, 4))
    plt.bar(categories, values, color=colors, alpha=alpha)
    
    if xlabel:
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)
    if title:
        plt.title(title)
    if ylim:
        plt.ylim(*ylim)
    if grid:
        plt.grid(axis="y")
    plt.tight_layout()
    if show:
        plt.show()