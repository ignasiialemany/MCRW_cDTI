# Monte Carlo Random Walk (MCRW) Simulation for Diffusion MRI

A C++ implementation of Monte Carlo simulations for diffusion MRI in biological tissues, with a focus on cardiac muscle tissue under different physiological states.

## Overview

This software simulates water diffusion in complex tissue environments using a Monte Carlo approach. It models:

- Water molecule movement in intracellular and extracellular spaces
- Membrane permeability effects on diffusion
- Tissue deformation and strain
- MRI pulse sequence effects

The simulation generates realistic diffusion MRI signals by tracking the phase accumulation of particles as they diffuse through the substrate during an MRI sequence.

## Quick Start with Docker

The simulation can be run easily using Docker, which handles all dependencies automatically. For some operations, you may need to use two terminals.

### Step 1: Build the Docker Image

Open your first terminal:

```bash
# Clone the repository (if you haven't already)
git clone https://github.com/ignasiialemany/MCRW_cDTI.git
cd MCRW_cDTI

# Build the Docker image
docker build -t mcrw .
```

### Step 2: Run the Container in Two Terminals

#### Terminal 1: Run the Simulation Container

```bash
# Run the container with mounted volume for simulation
docker run -it --rm -v $(pwd):/app mcrw
```

Inside this terminal, you'll compile and run the simulation:

```bash
# Create build directory
mkdir -p build/Release && cd build/Release

# Configure with CMake
cmake ../..

# Compile
make -j$(nproc)

# Run the simulation
./3DRandomwalk
```

## Using the SimulationRunner Class

The `SimulationRunner` class provides a simplified interface to run the simulation. The basic workflow is:

1. **Initialize**: Create a SimulationRunner instance with the number of particles, random seed, and strain function.
2. **Setup Sequence**: Load MRI sequence parameters from a YAML file.
3. **Setup Geometry**: Load tissue geometry from either MAT or OFF files.
4. **Run Simulation**: Configure simulation parameters and execute.

Example code:

```cpp
// Initialize SimulationRunner
int N_PARTICLES = 10000;
int SEED = 1;
std::function<double(double)> strain_function = [](double t) { return 0.0; };
auto runner = SimulationRunner(N_PARTICLES, SEED, strain_function);

// Setup sequence
runner.setupSequence("path/to/sequence.yaml");

// Setup geometry
double angle = 0.01;  // degrees per micrometer in axis Y
bool shift_block = false;
Eigen::VectorXd voxel(6);
voxel << min_x, max_x, min_y, max_y, min_z, max_z;
Eigen::VectorXd block_size(3);
block_size << dx, dy, dz;

// Option 1: MAT file geometry
runner.setupGeometryMatfile(angle, shift_block, "path/to/geometry.mat", voxel, block_size);

// Option 2: OFF file geometry
std::vector<std::string> geometry_paths = {"path/to/geometry1.off", "path/to/geometry2.off"};
runner.setupGeometryOffFiles(angle, shift_block, geometry_paths, voxel, block_size);

// Run simulation
parameters params;
params.isOutput = false;     // Save particle positions at each time step
params.kappa = 0.0;          // Membrane permeability
params.D_ecs = 2.5;          // Extracellular diffusivity (μm²/ms)
params.D_ics = 1.0;          // Intracellular diffusivity (μm²/ms)
params.cores = 64;           // Number of CPU cores to use
params.isDeformed = false;   // Apply strain deformation
params.strain_step_size = 100; // Time step for strain updates (ms)

runner.runSimulation(params, "output_filename");
```

## Output

The simulation produces CSV files containing:
- Initial particle positions (`output_filename_init.csv`)
- Final particle positions (`output_filename_final.csv`)

## Code Structure

The codebase is organized into several main components:

- **MonteCarlo/**: Core simulation engine and particle tracking
- **Substrate/**: Tissue environment definitions and geometry handling
- **Geometry/**: Geometric primitives and operations for tissue modeling
- **MRI/**: MRI sequence parameters and pulse definitions

## Dependencies

- Eigen (linear algebra)
- CGAL (computational geometry)
- Boost (filesystem, etc.)
- yaml-cpp (configuration parsing)
- OpenMP (parallelization)

All dependencies are automatically installed when using Docker.

## Advanced Configuration

For advanced users who need to modify the simulation settings beyond the parameters provided by SimulationRunner:

### Geometry Configuration

- **Block geometry**: Represents cardiac muscle cells (myocytes)
- **Voxel definition**: Represents the MRI imaging volume with buffer zones for diffusing particles
- **Strain function**: Models tissue deformation during cardiac cycles

### Simulation Parameters

Key parameters that can be adjusted include:
- `kappa`: Membrane permeability (0 = impermeable, higher values = more permeable)
- `D_ecs`: Extracellular diffusivity (μm²/ms)
- `D_ics`: Intracellular diffusivity (μm²/ms)
- `cores`: Number of CPU cores to use
- `isDeformed`: Whether to apply strain deformation
- `strain_step_size`: Time step for strain updates (ms)