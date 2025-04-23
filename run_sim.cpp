/**
 * @file run_sim.cpp
 * @brief Monte Carlo Random Walk (MCRW) simulation for diffusion MRI
 * 
 * This program simulates water diffusion in biological tissues using a Monte Carlo approach.
 * It models particles moving through different tissue environments (intracellular and extracellular)
 * with different diffusion properties and membrane permeability.
 * 
 * The program is part of a larger codebase that simulates diffusion MRI in biological tissues.
 * The main components include:
 * - MonteCarlo: Contains the simulation and walkers classes for particle movement
 * - Substrate: Defines the geometry and properties of the tissue environment
 * - Geometry: Handles the geometric primitives and operations
 * - MRI: Implements MRI sequence parameters and b-value calculations
 * 
 * This file orchestrates the entire simulation process:
 * 1. Loads tissue geometry from files
 * 2. Sets up transformation parameters for tissue deformation
 * 3. Configures the MRI sequence
 * 4. Defines strain functions for cardiac motion
 * 5. Initializes the substrate with appropriate geometry
 * 6. Seeds particles in the simulation space
 * 7. Executes the Monte Carlo simulation
 * 8. Writes output data to CSV files
 * 
 * Command line arguments:
 * argv[1]: Tissue state ("relaxed" or "compressed")
 * argv[2]: Extracellular volume fraction (ECV)
 * argv[3]: Whether to apply strain (0 for no, 1 for yes)
 * argv[4]: Index for permeability value (kappa) to use from predefined list
 * argv[5]: Random seed
 */

#include <iostream>
#include "yaml-cpp/yaml.h"
#include "MonteCarlo/simulation.h"
#include <boost/filesystem.hpp>
#include <cmath>
#include <fstream>

/**
 * @brief Writes two vectors to a CSV file
 * 
 * @param filename Output CSV filename
 * @param x First column data vector
 * @param y Second column data vector
 */
void writeCSVFile(const std::string& filename, const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
    std::ofstream file(filename);
    if (file.is_open()) {
        for (int i = 0; i < x.size(); ++i) {
            file << x(i) << "," << y(i) << std::endl;
        }
        file.close();
    } else {
        std::cerr << "Unable to open file " << filename << std::endl;
    }
}

/**
 * @brief Computes the cumulative sum of vector elements
 * 
 * @param x Input vector
 * @return Eigen::VectorXd Vector containing cumulative sums
 */
Eigen::VectorXd computeCumulativeArray(const Eigen::VectorXd& x) {
    Eigen::VectorXd cumulative(x.size());
    cumulative(0) = x(0);
    for (int i = 1; i < x.size(); ++i) {
        cumulative(i) = cumulative(i - 1) + x(i);
    }
    return cumulative;
}

/**
 * @brief Computes the cumulative integral of a function defined by x,y pairs
 * 
 * Uses the trapezoidal rule for numerical integration.
 * 
 * @param x X-coordinates (typically time points)
 * @param y Y-coordinates (function values)
 * @return Eigen::VectorXd Vector containing cumulative integral values
 */
Eigen::VectorXd cumulativeIntegral(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
    // Ensure that x and y are of the same size
    assert(x.size() == y.size());

    Eigen::VectorXd integral(y.size());
    integral(0) = 0; // Start with an integral of 0

    // Perform the cumulative integration using the trapezoidal rule
    for (int i = 1; i < x.size(); ++i) {
        double dx = x(i) - x(i - 1);
        double avgHeight = (y(i) + y(i - 1)) / 2.0;
        integral(i) = integral(i - 1) + dx * avgHeight;
    }

    return integral;
}

/**
 * @brief Computes the definite integral of a function defined by x,y pairs
 * 
 * Uses the trapezoidal rule for numerical integration.
 * 
 * @param x X-coordinates (typically time points)
 * @param y Y-coordinates (function values)
 * @return double Total integral value
 */
double computeIntegral(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
    // Ensure that x and y are of the same size
    assert(x.size() == y.size());

    double integral = 0.;
    // Perform the integration using the trapezoidal rule
    for (int i = 1; i < x.size(); ++i) {
        double dx = x(i) - x(i - 1);
        double avgHeight = (y(i) + y(i - 1)) / 2.0;
        integral = integral + dx * avgHeight;
    }

    return integral;
}

int main(int argc, char *argv[])
{
    // Process command line arguments
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " <state> <ecv> <strain_flag> <kappa_index> <seed>" << std::endl;
        return 1;
    }
    
    // Set up file paths
    auto current_path = boost::filesystem::current_path();
    auto parent_path = current_path.parent_path();
    std::string parent_path_str = parent_path.string();
    // Get substring of parent_path_str subtracting the last 5 characters
    std::string grandparent_path_str = parent_path_str.substr(0, parent_path_str.size() - 5);
    std::string file_path = "geometry_1.mat";
    std::string seq_path = "sequence.yaml";
    std::string full_path = grandparent_path_str + file_path;
    std::string seq_full_path = grandparent_path_str + seq_path;

    // Parse command line arguments
    std::string state = argv[1];  // Tissue state: "relaxed" or "compressed"
    std::string ecv = argv[2];    // Extracellular volume fraction
    bool isStrain = std::stoi(argv[3]);  // Whether to apply strain
    int index_kappa = std::stoi(argv[4]);  // Index for permeability value
    int seed = std::stoi(argv[5]);  // Random seed

    // Define available permeability (kappa) values in μm/ms
    std::vector<double> kappas = {0.,0.025,0.05, 0.1, 0.25, 0.5, 1., 2.5, 5., 10., 25.,50.,100.,1000., 100000.};
    
    // Load geometry data from MAT file
    auto data = utility_substrate::read_mat_file(full_path);

    // Configure transformation settings
    bool isIdentity = false;
    bool shift_block = true;
    double angle = 0.;

    // Initialize transformation with appropriate block dimensions based on tissue state
    transform t(angle, shift_block, isIdentity);
    if (state == "compressed"){
        t.set_block(31.65, 58.14, 91.84);
    }
    else{
        t.set_block(28.66, 52.66, 112);
    }
    
    // Load sequence parameters from YAML file
    sequence seq("../../sequence.yaml");
    seq.create();

    // Calculate cumulative time points and export sequence data
    Eigen::VectorXd cumulative_dt = computeCumulativeArray(seq.dt);
    std::cout << seq.dt << std::endl;
    std::cout << seq.gG << std::endl;
    writeCSVFile("dt_gG.csv", cumulative_dt, seq.gG);
    
    // Calculate b-value for the sequence
    auto cumulative_integral = cumulativeIntegral(cumulative_dt, seq.gG);
    Eigen::VectorXd squared_cumulative_integral = cumulative_integral.array().square();
    auto bvalue = computeIntegral(cumulative_dt, squared_cumulative_integral);
    std::cout << "bvalue: " << bvalue << std::endl;

    // Define strain functions for different physiological states
    
    // No strain - constant zero function
    std::function<double(double)> no_strain = [](double t){
        return 0.;
    };

    // End-systolic strain - linear ramp up and down
    std::function<double(double)> end_systolic_strain = [](double t){
        if (t < 500.){
            return (t/500.) * (0.2195);  // Linear increase to max strain
        }
        else{
            return 0.2195 - (t-500.)/500. * (0.2195);  // Linear decrease back to zero
        }
    };

    // End-diastolic strain - negative strain profile
    std::function<double(double)> end_diastolic_strain = [](double t){
        if (t < 500.){
            return t/500. * -0.18;  // Linear decrease to minimum strain
        }
        else if (t < 1000){
            return (0.18/500.) * (t-500.) - 0.18;  // Linear increase back to zero
        }
        else{
            return 0.;  // No strain after t=1000
        }
    };

    // Sinusoidal strain - sine wave pattern
    std::function<double(double)> sinusoidal_strain = [](double t){
        if (t <= 500){
            return -0.18 * std::sin(t/500. * M_PI);  // First half of sine wave
        }
        else if (t < 1000){
            return 0.18 * std::sin((t-500.)/500. * M_PI);  // Second half of sine wave
        }
        else{
            return 0.;  // No strain after t=1000
        }
    };

    // Load substrate geometry from OFF files
    std::vector<std::string> polygons;
    for (int i = 0; i <= 8; i++) {
        std::string path = "../../" + state + "/" + "simple_domain_ECV_" + ecv + "/" + "cuboid_" + std::to_string(i) + ".off";
        std::cout << path << std::endl;
        polygons.push_back(path);
    }
    
    // Select appropriate strain function based on input parameters
    int index_strain = 0;  // Default to no strain
    std::string strain_name = "no_strain";
    std::vector<std::function<double(double)>> strain_functions = {
        no_strain, end_systolic_strain, end_diastolic_strain, sinusoidal_strain
    };

    if (isStrain) {
        strain_name = "strained";
        if (state == "compressed") {
            index_strain = 1;  // End-systolic strain
        }
        else if (state == "relaxed") {
            index_strain = 2;  // End-diastolic strain
        }
        else if (state == "sinusoidal") {
            index_strain = 3;  // Sinusoidal strain
        }
    }

    // Initialize substrate with selected geometry and strain function
    substrate sub(polygons, t, strain_functions[index_strain]);

    // Define voxel dimensions with buffer zone
    Eigen::VectorXd voxel(6);
    double buffer_zone = std::sqrt(6 * 2.5 * seq.dt.sum());  // Calculate buffer based on diffusion distance
    voxel << 0.0 - buffer_zone, 2800 + buffer_zone, 
            0.0 - buffer_zone, 2800 + buffer_zone, 
            0.0 - buffer_zone, 2800 + buffer_zone;
    sub.setVoxel(voxel);

    // Create time steps for strain calculations
    double total_time = seq.dt.sum();
    double step = 100.;  // Time step size in ms
    int numElements = static_cast<int>(std::floor(total_time / step)) + 1;
    Eigen::VectorXd strain_array_time(numElements);
    
    // Fill time steps array, ensuring total matches sequence duration
    for (int i = 0; i < numElements - 1; ++i) {
        strain_array_time(i) = step;
    }
    strain_array_time(numElements - 1) = total_time - strain_array_time.head(numElements - 1).sum();
    std::cout << strain_array_time << std::endl;

    // Initialize simulation with particle count and random seed
    simulation sim(50000, seed);

    // Configure simulation parameters
    double kappa = kappas[index_kappa];
    sim.params.isOutput = false;
    sim.params.kappa = kappa;  // Membrane permeability (μm/ms)
    sim.params.D_ecs = 2.5;    // Extracellular diffusion coefficient (μm²/ms)
    sim.params.D_ics = 1.0;    // Intracellular diffusion coefficient (μm²/ms)
    sim.params.cores = 60;     // Number of CPU cores to use
    sim.params.isDeformed = isStrain;

    // Precompute substrate geometry for all time points
    auto start = std::chrono::high_resolution_clock::now();
    sim.precomputeSubstrate(sub, strain_array_time);
    auto end = std::chrono::high_resolution_clock::now();
    
    // Seed particles in the simulation volume
    sim.seedParticlesInBox(sub);
    // Alternative seeding method (commented out):
    // sim.seedingInMyocytes(sub);
    
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Precomputed substrate took: " << elapsed.count() << " s\n";
    std::cout << "cores: " << sim.params.cores << std::endl;

    // Start simulation and track initial particle positions
    std::cout << "Starting simulation" << std::endl;
    std::string file_init_positions = "kappa_" + std::to_string(kappa) + 
                                      "_state_" + state + 
                                      "_ecv_" + ecv + 
                                      "_" + strain_name + 
                                      "_init.csv";
    sim.writeParticlesState(file_init_positions, sub);

    // Run the simulation
    auto start_sim = std::chrono::high_resolution_clock::now();
    sim.performScan(sub, seq);
    auto end_sim = std::chrono::high_resolution_clock::now();   
    std::chrono::duration<double> elapsed_sim = end_sim - start_sim;
    std::cout << "Simulation took: " << elapsed_sim.count() << " s\n";

    // Write final particle positions to file
    std::string file_final_positions = "kappa_" + std::to_string(kappa) + 
                                       "_state_" + state + 
                                       "_ecv_" + ecv + 
                                       "_" + strain_name + 
                                       "_final.csv";
    sim.writeParticlesState(file_final_positions, sub);
    std::cout << "Finished simulation" << std::endl;
    
    return 0;
}

