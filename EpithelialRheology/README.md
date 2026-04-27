# EpithelialRheology — C++ Chaste User Project

Direct C++ port of the PyChaste experiment7.6 simulation, structured as a
reusable foundation for epithelial tissue rheology work.

---

## Project Structure

```
EpithelialRheology/
├── CMakeLists.txt                                   # Chaste build integration
├── src/
│   ├── NeighborDependentProliferationModifier.hpp   # Custom modifier (header)
│   └── NeighborDependentProliferationModifier.cpp   # Custom modifier (impl)
├── test/
│   ├── TestEpithelialFillSimulation.hpp             # Main test (port of 7.6)
│   └── ContinuousTestPack.txt                       # Test registration
└── README.md
```

## Setup

### 1. Install Chaste (if not already)

Follow the official guide:
https://chaste.github.io/docs/installguides/

The recommended route is the **Docker image** or **building from source** on
Ubuntu. You need the `cell_based` component.

### 2. Add this project to Chaste

Copy or symlink the `EpithelialRheology/` folder into Chaste's `projects/`
directory:

```bash
# From the Chaste source root:
ln -s /path/to/EpithelialRheology projects/EpithelialRheology
```

### 3. Build

```bash
cd <chaste-build-directory>
cmake <chaste-source-directory>
make -j$(nproc) project_EpithelialRheology
```

### 4. Run

```bash
ctest -R TestEpithelialFillSimulation -V
```

Output VTU files will be written to the Chaste test output directory
(typically `~/testoutput/EpithelialRheology/FillSimulation/`).

---

## What's Ported from experiment7.6.py

| Feature                            | Python (PyChaste)                     | C++ (this project)                              |
|------------------------------------|---------------------------------------|-------------------------------------------------|
| Vertex mesh (honeycomb)            | `HoneycombVertexMeshGenerator`        | Same class, native API                          |
| Contact inhibition cell cycle      | `ContactInhibitionCellCycleModel`     | Same class, full parameter control              |
| Neighbor-based proliferation block | Manual loop before each `Solve()`     | `NeighborDependentProliferationModifier` (runs every dt) |
| NagaiHonda force                   | Pre-wrapped binding                   | Native class, same parameters                   |
| Plane boundary conditions          | 4× `PlaneBoundaryCondition`           | Same, native API                                |
| Coverage-based early stopping      | Incremental `Solve()` + area check    | Same incremental pattern                        |
| T3 swap disable + threshold tuning | `SetCheckForT3Swaps(False)`           | Same                                            |
| Fixed RNG seed                     | `RandomNumberGenerator.Reseed()`      | Same                                            |

### Key improvement over PyChaste

The `NeighborDependentProliferationModifier` is a proper Chaste simulation
modifier with:
- Serialization support (checkpointing works)
- Runs at every timestep (not just once per hour like the Python loop)
- Configurable threshold parameter
- Reusable across different test configurations

---

## Extending for Rheology

This is your starting point. Here's what you can now build **that PyChaste
couldn't give you**:

### Custom forces (stress extraction)
Create a subclass of `AbstractForce<2>` to compute and log per-cell stress
tensors. Put the header/source in `src/` — the build system picks them up
automatically.

```cpp
// src/StressTrackingForce.hpp — skeleton
#include "AbstractForce.hpp"
template<unsigned DIM>
class StressTrackingForce : public AbstractForce<DIM>
{
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation) override;
    // ... compute & store stress per cell using CellData
};
```

### Custom cell cycle models
Subclass `AbstractCellCycleModel` to couple mechanics (cell area, stress)
directly to proliferation decisions — not possible in PyChaste.

### Proper apoptosis
Subclass `AbstractCellKiller<2>` for mechanically-triggered cell death.
No more try-except workarounds.

### Curved boundary conditions
Subclass `AbstractCellPopulationBoundaryCondition<2,2>` for native circle
or ellipse constraints — no more 10-plane approximations.

### Custom output writers
Subclass `AbstractCellWriter` or `AbstractCellPopulationWriter` to output
rheological quantities (strain rate, local viscosity, etc.) at each timestep.

---

## Adding New Tests

1. Create a new `.hpp` file in `test/` (e.g., `TestRheologyMeasurement.hpp`)
2. Add the filename to `test/ContinuousTestPack.txt`
3. Rebuild: `make -j$(nproc) project_EpithelialRheology`
4. Run: `ctest -R TestRheologyMeasurement -V`
