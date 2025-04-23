# Monte Carlo Random Walk (MCRW) Simulation for Diffusion MRI

A C++ implementation of Monte Carlo simulations for diffusion MRI in biological tissues, with a focus on cardiac muscle tissue under different physiological states.

## Overview

This software simulates water diffusion in complex tissue environments using a Monte Carlo approach. It models:

- Water molecule movement in intracellular and extracellular spaces
- Membrane permeability effects on diffusion
- Tissue deformation and strain
- MRI pulse sequence effects

The simulation generates realistic diffusion MRI signals by tracking the phase accumulation of particles as they diffuse through the substrate during an MRI sequence.

## Code Structure

The codebase is organized into several main components:

- **MonteCarlo/**: Core simulation classes
  - `simulation.h/cpp`: Main simulation engine
  - `walkers.h/cpp`: Particle management and tracking
  
- **Substrate/**: Tissue environment definitions
  - `substrate.h/cpp`: Substrate geometry and properties
  - `transform.h/cpp`: Transformation and strain handling
  - `utility_substrate.h/cpp`: Utility functions for substrate operations
  
- **Geometry/**: Geometric primitives and operations
  - `polyhedronSet.h/cpp`: Collection of 3D shapes for strain transformation
  
- **MRI/**: MRI sequence related code
  - `sequence.h/cpp`: Pulse sequence definitions

- **run_sim.cpp**: Main entry point that orchestrates the simulation

## Building and Requirements

### Dependencies

- Eigen (linear algebra)
- CGAL (computational geometry)
- Boost (filesystem, etc.)
- yaml-cpp (configuration parsing)
- OpenMP (parallelization)

### Build Instructions

```bash
mkdir build
cd build
cmake ..
make
```

## Usage

The simulation accepts the following command-line arguments:

```bash
./run_sim <state> <ecv> <strain_flag> <kappa_index> <seed>
```

Where:
- `state`: Tissue state ("relaxed" or "compressed")
- `ecv`: Extracellular volume fraction
- `strain_flag`: Whether to apply strain (0 for no, 1 for yes)
- `kappa_index`: Index for permeability value (kappa) from the predefined list
- `seed`: Random seed for the simulation

## Example

```bash
./run_sim relaxed 0.2 1 5 42
```

This runs a simulation with:
- Relaxed tissue state
- 20% extracellular volume
- With strain applied
- Permeability index 5 (corresponding to kappa = 1.0)
- Random seed 42

## Output

The simulation produces CSV files containing:
- Initial particle positions
- Final particle positions
- Diffusion signal data

Output filenames include the simulation parameters for easy identification.

## Geometric Representation

The simulation uses a sophisticated geometry system to model tissue microstructure and particle movement:

### Block Geometry

The block geometry represents the cardiac muscle cells (myocytes):

- **Myocyte Representation**:
  - Each myocyte is modeled as a rectangular cuboid with dimensions (dx, dy, dz)
  - Different dimensions are used based on physiological state:
    - Compressed: (31.65, 58.14, 91.84) μm
    - Relaxed: (28.66, 52.66, 112) μm

- **polyhedronSet Class**:
  - Manages collections of 3D shapes that form the tissue structure
  - Stores the geometry at different time points during deformation
  - Based on CGAL triangulations for efficient geometric operations

- **Tissue Geometry Loading**:
  - OFF files (Object File Format) define the 3D tissue geometry
  - Multiple cuboids are combined to create the complete tissue structure
  - Different geometries are provided for various extracellular volume fractions (ECV)

### Voxel Definition

The voxel represents the MRI imaging volume:

- **Structure**:
  - Defined by 6 coordinates: [minX, maxX, minY, maxY, minZ, maxZ]
  - Default size covers 0 to 2800 μm in each dimension

- **Buffer Zone**:
  - The voxel includes a buffer zone to capture diffusing particles
  - Calculated as sqrt(6 × D × T) where:
    - D is the diffusion coefficient (2.5 μm²/ms)
    - T is the total simulation time
  - Ensures accurate representation of diffusion effects at boundaries

### Coordinate Systems

The simulation uses multiple coordinate systems for accurate modeling:

- **Global Coordinates**: Position in the entire simulation space
- **Local Coordinates**: Position relative to a specific myocyte
- **Deformed Coordinates**: Positions accounting for tissue strain

The `transform` class handles conversions between these coordinate systems and applies strain-induced deformations to the geometry.

### Particle-Geometry Interaction

- **Compartment Detection**:
  - The simulation tracks whether particles are in intracellular or extracellular space
  - Uses spatial indexing for efficient compartment determination

- **Membrane Interaction**:
  - When particles encounter cell membranes, permeability (kappa) determines crossing probability
  - Permeability values range from 0 (impermeable) to 100,000 μm/s (highly permeable)
  - Particles may reflect off membranes or cross them based on probabilistic rules

- **Strain Effects**:
  - For cardiac motion simulations, geometry deforms over time
  - Precomputed geometric configurations are used at different time points
  - Strain functions model different cardiac states:
    - End-systolic strain (compression)
    - End-diastolic strain (relaxation)
    - Sinusoidal strain (cyclic contraction)
