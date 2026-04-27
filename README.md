# BIOE242 Final Project — Epithelial Tissue Mechanics via Graph Neural Networks

**Christian Chung | BIOE242 Spring 2026**

## Overview

Can cell-level contact topology encode the mechanical state of an epithelial tissue? This project addresses that question by generating synthetic epithelial tissues via Chaste vertex model simulations, extracting cell contact graphs, and training graph neural networks (GNNs) to classify tissue mechanical state and perturbed biological parameters.

Two classification tasks are explored:
- **Binary:** fluid vs. solid phase (shape index threshold q* = 3.81)
- **3-class suite:** identify which biological parameter was perturbed — Adhesion, Deformation Energy, or Proliferation

---

## Repository Structure

```
bioe242-final-project-cc/
├── BIOE242_Checkpoint1_DataPreparation.ipynb           # Data prep — single seed
├── BIOE242_Checkpoint1_DataPreparation_Replicated.ipynb # Data prep — 3 seeds
├── BIOE242_Checkpoint2_ModelConstruction.ipynb         # Baseline models & unsupervised analysis
├── BIOE242_Checkpoint3_Regularization.ipynb            # Hyperparameter tuning
├── BIOE242_Checkpoint4_ProductionMode.ipynb            # Production model & temporal trajectory
├── BIOE242_Checkpoint5_FinalResults.ipynb              # Final analysis & conclusions
│
├── BIOE242 Data Bank/                                  # Processed data (single seed, 2,793 graphs)
├── BIOE242 Data Bank Complete/                         # Processed data (3 seeds, 8,605 graphs)
│
├── EpithelialRheology/                                 # C++ Chaste user project (simulation source)
│   ├── src/
│   │   ├── NeighborDependentProliferationModifier.hpp
│   │   └── NeighborDependentProliferationModifier.cpp
│   └── CMakeLists.txt
│
├── replicated_testoutput/                              # Raw Chaste simulation output (.vtu files)
└── improvement_journey.png                             # Training accuracy improvement visualization
```

---

## Methods

### 1. Simulation (Chaste Vertex Model)

Epithelial tissue dynamics were simulated using the [Chaste](https://chaste.cs.ox.ac.uk/) vertex model framework with the FarhadifarForce (line tension + area elasticity). A 15×4 honeycomb lattice was used. Three biological parameters were swept independently across 10 values each:

| Suite | Parameter | Biological interpretation |
|---|---|---|
| Adhesion | Line tension coefficient | Cell-cell adhesion strength |
| DeformationEnergy | Area elasticity coefficient | Resistance to cell shape change |
| Proliferation | Cell cycle scale | Rate of cell division |

Three independent seeds were run per condition for a total of 93 simulation runs. A custom `NeighborDependentProliferationModifier` was implemented to arrest interior cells with ≥6 neighbors, preventing overcrowding.

### 2. Graph Construction (Checkpoint 1)

Each simulation snapshot is parsed from `.vtu` mesh files using `meshio`. A cell contact graph is built where:
- **Nodes** = individual cells with features: `[area, perimeter, shape_index, n_neighbors]`
- **Edges** = shared cell boundaries (undirected)

To prevent temporal autocorrelation and class imbalance, each simulation variant is capped at N_MAX = 300 snapshots. Train/val/test split is done at the variant level (7/2/1) to prevent data leakage. Node features are Z-score normalized using training set statistics.

**Final dataset size (replicated):** 8,605 graphs — 5,984 train / 1,737 val / 884 test

### 3. Models (Checkpoint 2)

Two model families are benchmarked:

- **MLP (topology-blind null baseline):** operates on 8-dimensional pooled features (mean + std of node features across all cells)
- **GraphSAGE GNN:** 3-layer message-passing network with mean aggregation, operates on cell contact graph

Both are evaluated on binary (fluid/solid) and 3-class (suite) classification.

### 4. Hyperparameter Tuning (Checkpoint 3)

Sequential sweep over five hyperparameters for the GNN:

| Hyperparameter | Range | Best value |
|---|---|---|
| Dropout | 0.0 – 0.5 | 0.0 |
| Weight decay | 0.0 – 5e-3 | 1e-3 |
| GNN depth (layers) | 2 – 5 | 5 |
| Hidden dimension | 64 – 256 | 128 |
| Learning rate | 5e-4 – 2e-3 | 5e-4 |

A critical methodological fix — **best-checkpoint saving** (save model whenever val accuracy improves, restore before evaluation) — eliminated severe overfitting that was masking model quality in baseline experiments.

---

## Results

### Binary Classification (Fluid vs. Solid)

Both the MLP and GNN achieve >99% test accuracy. Shape index (q) is a near-perfect scalar order parameter, consistent with vertex model jamming theory (q* = 3.81).

### 3-Class Suite Classification

| Model | Test Accuracy |
|---|---|
| Baseline GNN | 37.8% |
| Baseline MLP | 45.7% |
| Tuned GNN | **60.4%** |
| Tuned MLP | 50.2% |

The tuned GNN's +10 pp improvement over the tuned MLP demonstrates that cell contact topology carries genuine mechanistic signal beyond bulk statistics. The ~60% ceiling reflects real topological overlap across suites, especially for the Proliferation condition.

**Per-class performance (tuned GNN):**
- Adhesion: **84.5%** — line tension creates distinctive polygon distributions
- DeformationEnergy: **56.3%** — moderate cell shape differences
- Proliferation: **41.7%** — transient topology disruption from cell division makes this hardest

### Temporal Trajectory (Checkpoint 4)

The production GNN (retrained on all 8,305 graphs) was applied to 209 consecutive snapshots of the fastest proliferating simulation. Predicted class probabilities track the correct class throughout the simulation and correlate with tissue-mean shape index evolution, validating the model's biological interpretability.

---

## Dependencies

No `requirements.txt` is included. The project uses:

- **Python:** PyTorch, PyTorch Geometric, NumPy, Pandas, scikit-learn, matplotlib, meshio, UMAP
- **C++ simulation:** [Chaste](https://chaste.cs.ox.ac.uk/) (for running new simulations; not required to run notebooks on existing data)

---

## Limitations & Future Work

- No external validation (all data is synthetic Chaste output)
- Class imbalance not addressed (weighted loss or oversampling not explored)
- CPU-only training limited hyperparameter sweep scope
- Single-snapshot analysis — temporal graph networks could capture dynamics directly
- More expressive GNNs (e.g., GIN) and edge features (edge length, tension) could close the MLP gap further
