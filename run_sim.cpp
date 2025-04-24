#include <iostream>
#include "SimulationRunner.h"
#include <boost/filesystem.hpp>

int main(int argc, char* argv[]) {

    //GLOBAL PARAMETERS
    int SEED = 1;
    int N_PARTICLES = 10000;

    std::function<double(double)> strain_function = [](double t) { return 0.0; };
    auto runner = SimulationRunner(N_PARTICLES, SEED, strain_function);

    boost::filesystem::path executable_path = boost::filesystem::current_path();
    boost::filesystem::path root_dir = executable_path.parent_path().parent_path();

    std::cout << "Root directory: " << root_dir << std::endl;
    std::cout << "Executable path: " << executable_path << std::endl;

    std::string sequence_path = root_dir.string() + "/sequence.yaml";
    std::string geometry_path = root_dir.string() + "/geometry_1.mat";

    //Setup sequence
    std::cout << "sequence path: " << sequence_path << std::endl;
    runner.setupSequence(sequence_path);
   
    //Geometry parameters
    double angle = 0.01; //degrees per micrometer in axis Y - default is 0.01 (10 degrees per mm)
    bool shift_block = false;
    //We assume D=3 (diffusivity coefficient of water at 37C) for the buffer size
    double buffer_size = std::sqrt(2*3*runner.sequence_ptr_->dt.sum());

    Eigen::VectorXd voxel(6);
    // X=8000, Y=2800, Z=2800 dimensions for cDTI voxel (Note that X is the longest dimension) in um
    voxel << 0 - buffer_size, 8000 + buffer_size, 0 - buffer_size, 2800 + buffer_size, 0 - buffer_size, 2800 + buffer_size;

    Eigen::VectorXd block_size(3);
    // delta_x=495.3992, delta_y=392.3432, delta_z=126.5612 dimensions for the block in um 
    block_size << 495.4, 392.4, 126.6;

    //Setup geometry from .mat file - we can also set up with .off files (see below)
    runner.setupGeometryMatfile(angle, shift_block, geometry_path, voxel, block_size);

    /*
    //Relaxed geometry
    std::vector<std::string> geometry_paths;
    geometry_paths.push_back(root_dir.string() + "/relaxed/geometry_2.off");
    geometry_paths.push_back(root_dir.string() + "/relaxed/geometry_3.off");
    geometry_paths.push_back(root_dir.string() + "/relaxed/geometry_4.off");

    //Setup geometry from .off files
    runner.setupGeometryOffFiles(angle, shift_block, geometry_paths, voxel, block_size);
    */

    //Run simulation
    parameters params;
    //Set to true if you want to save the particle positions in each time step
    params.isOutput = false;
    //Kappa parameter for the simulation
    params.kappa = 0.;
    //Diffusivity of the ECS
    params.D_ecs = 2.5;
    //Diffusivity of the ICS
    params.D_ics = 1;
    //Number of cores for the simulation
    params.cores = 64;
    //Set to false if there is no strain (faster no need to store deformed geometries)
    params.isDeformed = false;
    //Strain step size in ms
    params.strain_step_size = 100;

    //Run simulation
    runner.runSimulation(params, "test_simulation");
    
    return 0;
}
