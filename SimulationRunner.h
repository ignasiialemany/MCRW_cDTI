// SimulationRunner.h
#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <chrono>
#include <Eigen/Dense>
#include "MonteCarlo/simulation.h"
#include "yaml-cpp/yaml.h"
#include <boost/filesystem.hpp>
#include <cmath>
#include <fstream>
#include <boost/optional.hpp>
#include <memory>

// Forward declarations
class simulation;
class transform;
class substrate;
class sequence;
struct parameters;

class SimulationRunner {
public:
    // Constructor
    SimulationRunner(int N_particles, int seed, std::function<double(double)> strain_function);

    // Setup methods
    void setupGeometryMatfile(double angle, bool shift_block, std::string geometry_path, Eigen::VectorXd voxel, Eigen::VectorXd block_size);
    void setupGeometryOffFiles(double angle, bool shift_block, std::vector<std::string> geometry_paths, Eigen::VectorXd voxel, Eigen::VectorXd block_size);
    void setupSequence(std::string seq_path);

    // Simulation methods
    Eigen::VectorXd createStrainTimeArray(double step_size);
    void runSimulation(parameters params, std::string file_name);

    // Static utility functions
    static Eigen::VectorXd computeCumulativeArray(const Eigen::VectorXd& x);
    static Eigen::VectorXd cumulativeIntegral(const Eigen::VectorXd& x, const Eigen::VectorXd& y);
    static double computeIntegral(const Eigen::VectorXd& x, const Eigen::VectorXd& y);

    std::unique_ptr<sequence> sequence_ptr_;

private:
    // Member variables
    int N_particles_;
    int seed_;
    std::function<double(double)> strain_function_;

    // Unique pointers to class instances
    std::unique_ptr<simulation> simulation_ptr_;
    std::unique_ptr<transform> transform_ptr_;
    std::unique_ptr<substrate> substrate_ptr_;
    
};