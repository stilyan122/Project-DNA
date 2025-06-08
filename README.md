# Project DNA

Simulating the stochastic evolution of DNA under the Jukes–Cantor model and quantifying convergence of nucleotide frequencies to uniformity.

---

## 📖 Table of Contents

- [Introduction](#introduction)  
- [Features](#features)  
- [Repository Structure](#repository-structure)  
- [Installation](#installation)  
- [Usage](#usage)  
  - [Population-Level Simulation](#population-level-simulation)  
  - [Sequence-Level Sampling](#sequence-level-sampling)  
- [Analysis & Visualization](#analysis--visualization)  
- [Key Findings](#key-findings)  
- [Contributing](#contributing)  
- [License](#license)  

---

## Introduction

Under neutral evolution, each nucleotide (A, C, G, T) mutates to any other with equal probability according to the Jukes–Cantor (JC) model. **Project DNA** implements both population-level and sequence-level simulations to empirically measure how many generations are needed for base frequencies to approach uniformity, and applies bootstrap methods to quantify stochastic variability.

---

## Features

- ✅ **JC Markov-Chain Model**: Generate the transition matrix for any mutation rate μ.  
- ✅ **Sequence-Level Simulation**: Mutate each base at rate μ per generation.  
- ✅ **Distance Metrics**: Compute L₁ distance between observed base frequencies and uniform distribution.  
- ✅ **Bootstrap Analysis**: Estimate confidence intervals for “time-to-uniformity” over many replicates.  
- ✅ **Visualization Helpers**: Quick plotting functions for convergence curves and frequency histograms.  

---

## Repository Structure

```text
Project DNA/
├─ .git/                    
├─ README.md               ← This file  
└─ notebooks_code/
   ├─ bio_intro.py         ← Nucleotide basics & context  
   ├─ bio_structures.py    ← DNASequence and data structures  
   ├─ mutation_random.py   ← IID mutation demo & die analogy  
   ├─ markov_mutations.py  ← JC transition matrix & chain updates  
   ├─ hypothesis_testing.py← Bootstrap & statistical tests  
   ├─ plotting.py          ← Quick plotting utilities  
   └─ Project_DNA.ipynb     ← Master notebook: runs end-to-end  
```

---

## Installation

1. **Clone the repo**  
   ```bash
   git clone https://github.com/<your-username>/project-dna.git
   cd project-dna
   ```
2. **Create & activate a virtual environment**  
   ```bash
   python3 -m venv venv
   source venv/bin/activate
   ```
3. **Install dependencies**  
   ```bash
   pip install -r requirements.txt
   ```
   *Requirements include:*  
   - Python 3.10+  
   - NumPy  
   - SciPy  
   - Matplotlib  
   - (Optional) NetworkX  

---

## Usage

### Population-Level Simulation

```python
from markov_mutations import jukes_cantor_matrix, is_stochastic_matrix
import numpy as np

μ = 0.01
P = jukes_cantor_matrix(mu=μ)
assert is_stochastic_matrix(P)

# Initial frequencies: all As
dist = np.array([1.0, 0.0, 0.0, 0.0])
for gen in range(1, 5001):
    dist = dist @ P
    # compute L1 distance to uniform [0.25, 0.25, 0.25, 0.25]
```

### Sequence-Level Sampling

```python
from bio_structures import DNASequence

seq = DNASequence.random(length=1_000_000, freq=[1,0,0,0])
for gen in range(1, 5001):
    seq.mutate(mu=μ)
    freqs = seq.base_frequencies()
    # compute L1 distance to uniform [0.25, 0.25, 0.25, 0.25]
```

---

## Analysis & Visualization

- **Convergence Curves**: Use `plotting.quick_line_plot()` to display L₁ distance over generations.  
- **Frequency Histograms**: Use `plotting.quick_bar_plot()` at snapshots (e.g. gen 0, gen 1000, gen 3000).  
- **Bootstrap**: Run `hypothesis_testing.bootstrap_time_to_uniform()` to get confidence intervals for threshold crossing.

---

## Key Findings

- **Convergence Speed**  
  For μ ≈ 0.01, the L₁ distance falls below 0.05 around **2,500–3,000 generations**.  
- **Stochastic Variability**  
  Bootstrap replicates show a spread of ±200–400 generations in “time-to-uniformity.”

---

## Contributing

1. Fork the repo  
2. Create a feature branch (`git checkout -b feature/XYZ`)  
3. Commit your changes (`git commit -m "Add XYZ"`)  
4. Push to your branch (`git push origin feature/XYZ`)  
5. Open a Pull Request

Please follow the existing code style, add tests where appropriate, and update this README if you add new functionality.

---

## License

This project is licensed under the [MIT License](LICENSE). Feel free to use, modify, and distribute under the terms of MIT.

---
