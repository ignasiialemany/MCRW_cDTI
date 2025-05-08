# Monte Carlo Random Walk (MCRW) Simulation for Diffusion MRI

A C++ implementation of Monte Carlo simulations for diffusion MRI in biological tissues, with a focus on cardiac muscle tissue under different physiological states.

## Overview

This software simulates water diffusion in complex tissue environments using a Monte Carlo approach. It models:

- Water molecule movement in intracellular and extracellular spaces
- Membrane permeability effects on diffusion
- Tissue deformation and strain
- MRI pulse sequence effects

The simulation generates realistic diffusion MRI signals by tracking the phase accumulation of particles as they diffuse through the substrate during an MRI sequence.

## Running the Simulation

There are three different ways to run this simulation code, depending on your needs:

### 1. Local Installation (for Development)

If you're actively developing the code, you can install all dependencies locally:

1. **Install Dependencies**: 
   - Install Anaconda or Miniconda first
   - Required packages: Eigen3, CGAL, Boost, yaml-cpp, OpenMP, Ceres, matioCpp, Catch2

2. **Configure CMake**:
   - Open `CMakeLists.txt` and uncomment the Anaconda environment configuration section:
   ```cmake
   # Uncomment these lines:
   set(CONDA_ENV_NAME "randomwalk")
   set(CONDA_PREFIX "~/anaconda3/envs/${CONDA_ENV_NAME}")
   set(CONDA_INCLUDE_DIR "${CONDA_PREFIX}/include")
   set(CONDA_LIBRARY_DIR "${CONDA_PREFIX}/lib")
   set(CMAKE_PREFIX_PATH ${CONDA_PREFIX} ${CMAKE_PREFIX_PATH})
   ```

3. **Build and Run**:
   ```bash
   mkdir -p build/Release && cd build/Release
   cmake ../..
   make -j$(nproc)
   ./3DRandomwalk
   ```

> **Note**: This setup process can be complex and may require adjustments based on your operating system and local environment.

### 2. Docker Container (Recommended for Most Users)

Docker handles all dependencies automatically, making it the easiest way to run the simulation:

```bash
# Clone the repository (if you haven't already)
git clone https://github.com/ignasiialemany/MCRW_cDTI.git
cd MCRW_cDTI

# Build the Docker image
docker build -t mcrw .

# Run the container with mounted volume
docker run -it -v $(pwd):/app mcrw
```

Once inside the container, you can run the simulation:

```bash
# Navigate to the executable
cd /app/build/Release

# Run the simulation with parameters
./3DRandomwalk --kappa 0.01 --strain-type diastolic
```

All your files in the host directory are synchronized with the `/app` directory in the container, so simulation outputs will be available on your local machine.

### 3. High-Performance Computing (HPC)

For running simulations on an HPC cluster:

1. **Setup Directory Structure**:
   Create and navigate to a directory structure like:
   ```
   root/
   ├── build/
   │   └── Release/
   ├── sequence.yaml
   └── geometry_1.mat
   ```

2. **Pull Singularity Image**:
   Inside the `root/build/Release` directory:
   ```bash
   singularity pull simulation.sif docker://ignasiialemany/mcrw-cdti:latest
   ```

3. **Run the HPC Script**:
   Copy the `scripts/hpc_example.sh` script to your HPC system and submit it using your cluster's job scheduler:
   ```bash
   qsub hpc_example.sh
   ```

The script will run the simulation using the Singularity container, with output files generated in the `root/build/Release` directory.

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

## Command Line Arguments

The simulation supports various command line arguments:

```
Options:
  --kappa VALUE             Set kappa value (default: 0.0)
  --seed VALUE              Set random seed (default: 1)
  --particles N             Set number of particles (default: 10000)
  --strain-type TYPE        Set strain type: diastolic, systolic, no_strain (default: diastolic)
  --angle VALUE             Set angle in degrees per micrometer (default: 0.01)
  --shift-block BOOL        Enable/disable shift block: true, false (default: false)
  --Decs VALUE              Set diffusivity of ECS (default: 2.5)
  --Dics VALUE              Set diffusivity of ICS (default: 1.0)
  --cores N                 Set number of cores (default: 64)
  --deformed BOOL           Enable/disable deformation: true, false (default: true)
  --strain-step-size VALUE  Set strain step size in ms (default: 100)
  --job-id ID               Job array ID for output naming
  --sequence-file PATH      Path to sequence file (default: [root_dir]/sequence.yaml)
  --geometry-file PATH      Path to geometry file (default: [root_dir]/geometry_1.mat)
```

## Output

The simulation produces CSV files containing:
- Initial particle positions (`output_filename_init.csv`)
- Final particle positions (`output_filename_final.csv`)

## Dependencies

- Eigen (linear algebra)
- CGAL (computational geometry)
- Boost (filesystem, etc.)
- yaml-cpp (configuration parsing)
- OpenMP (parallelization)
- matioCpp (MAT-file handling)
- Ceres (optimization)

All dependencies are automatically installed when using Docker or Singularity.