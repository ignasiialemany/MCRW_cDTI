#include <iostream>
#include <string>
#include <sstream>
#include "SimulationRunner.h"
#include <boost/filesystem.hpp>
#include <vector>

void printUsage(const char* programName) {
    std::cerr << "Usage: " << programName << " [options]\n"
              << "Options:\n"
              << "  --kappa VALUE             Set kappa value (default: 0.0)\n"
              << "  --seed VALUE              Set random seed (default: 1)\n"
              << "  --particles N             Set number of particles (default: 10000)\n"
              << "  --strain-type TYPE        Set strain type: diastolic, systolic, no_strain (default: diastolic)\n"
              << "  --angle VALUE             Set angle in degrees per micrometer (default: 0.01)\n"
              << "  --shift-block BOOL        Enable/disable shift block: true, false (default: false)\n"
              << "  --Decs VALUE              Set diffusivity of ECS (default: 2.5)\n"
              << "  --Dics VALUE              Set diffusivity of ICS (default: 1.0)\n"
              << "  --cores N                 Set number of cores (default: 64)\n"
              << "  --deformed BOOL           Enable/disable deformation: true, false (default: true)\n"
              << "  --strain-step-size VALUE  Set strain step size in ms (default: 100)\n"
              << "  --job-id ID               Job array ID for output naming\n"
              << "  --sequence-file PATH      Path to sequence file (default: [root_dir]/sequence.yaml)\n"
              << "  --geometry-file PATH      Path to geometry file (default: [root_dir]/geometry_1.mat)\n"
              << "  --help                    Display this help message\n";
}

// Helper function to convert string to bool
bool strToBool(const std::string& str) {
    if (str == "true" || str == "TRUE" || str == "1" || str == "yes" || str == "YES" || str == "y" || str == "Y") {
        return true;
    }
    return false;
}

int main(int argc, char* argv[]) {
    // Default parameter values
    double kappa = 0.0;
    int SEED = 1;
    int N_PARTICLES = 10000;
    std::string strain_type = "diastolic";
    double angle = 0.01;
    bool shift_block = true;
    double D_ecs = 2.5;
    double D_ics = 1.0;
    int cores = 64;
    bool isDeformed = true;
    double strain_step_size = 100.0;
    std::string job_id = "0";
    std::string sequence_file = "";  // Will be set based on root_dir if empty
    std::string geometry_file = "";  // Will be set based on root_dir if empty
    
    // Parse command-line arguments
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "--kappa" && i + 1 < argc) {
            kappa = std::stod(argv[++i]);
        }
        else if (arg == "--seed" && i + 1 < argc) {
            SEED = std::stoi(argv[++i]);
        }
        else if (arg == "--particles" && i + 1 < argc) {
            N_PARTICLES = std::stoi(argv[++i]);
        }
        else if (arg == "--strain-type" && i + 1 < argc) {
            strain_type = argv[++i];
            if (strain_type != "diastolic" && strain_type != "systolic" && strain_type != "no_strain") {
                std::cerr << "Invalid strain type. Must be 'diastolic', 'systolic', or 'no_strain'." << std::endl;
                return 1;
            }
        }
        else if (arg == "--angle" && i + 1 < argc) {
            angle = std::stod(argv[++i]);
        }
        else if (arg == "--shift-block" && i + 1 < argc) {
            shift_block = strToBool(argv[++i]);
        }
        else if (arg == "--Decs" && i + 1 < argc) {
            D_ecs = std::stod(argv[++i]);
        }
        else if (arg == "--Dics" && i + 1 < argc) {
            D_ics = std::stod(argv[++i]);
        }
        else if (arg == "--cores" && i + 1 < argc) {
            cores = std::stoi(argv[++i]);
        }
        else if (arg == "--deformed" && i + 1 < argc) {
            isDeformed = strToBool(argv[++i]);
        }
        else if (arg == "--strain-step-size" && i + 1 < argc) {
            strain_step_size = std::stod(argv[++i]);
        }
        else if (arg == "--job-id" && i + 1 < argc) {
            job_id = argv[++i];
        }
        else if (arg == "--sequence-file" && i + 1 < argc) {
            sequence_file = argv[++i];
        }
        else if (arg == "--geometry-file" && i + 1 < argc) {
            geometry_file = argv[++i];
        }
        else if (arg == "--help") {
            printUsage(argv[0]);
            return 0;
        }
        else {
            std::cerr << "Unknown argument: " << arg << std::endl;
            printUsage(argv[0]);
            return 1;
        }
    }
    
    std::cout << "Running simulation with parameters:" << std::endl;
    std::cout << "  Kappa: " << kappa << std::endl;
    std::cout << "  Seed: " << SEED << std::endl;
    std::cout << "  Particles: " << N_PARTICLES << std::endl;
    std::cout << "  Strain Type: " << strain_type << std::endl;
    std::cout << "  Angle: " << angle << std::endl;
    std::cout << "  Shift Block: " << (shift_block ? "true" : "false") << std::endl;
    std::cout << "  D_ecs: " << D_ecs << std::endl;
    std::cout << "  D_ics: " << D_ics << std::endl;
    std::cout << "  Cores: " << cores << std::endl;
    std::cout << "  Deformed: " << (isDeformed ? "true" : "false") << std::endl;
    std::cout << "  Strain Step Size: " << strain_step_size << std::endl;
    std::cout << "  Job ID: " << job_id << std::endl;

    // Select strain function based on strain_type parameter
    std::function<double(double)> strain_function;
    
    if (strain_type == "diastolic") {
        strain_function = [](double t){
            double result;
            if (t <= 300) {
                result = -0.15 * sin(2 * M_PI * t * 1 / 1200);
            } else if (t > 300 && t <= 800) {
                result = -0.15 * sin(2 * M_PI * (t + 200) * 1 / 2000);
            } else {
                result = 0.;
            }
            
            double scale = 10000; // Scale for 4 decimal places
            return std::round(result * scale) / scale;
        };
    } 
    else if (strain_type == "systolic") {
        strain_function = [](double t){
            double result;
            if (t < 500) {
                result = (1/0.85-1) * (1 - sin(2*M_PI*(t+500)*1/2000));
            }
            else if (t>=500 and t<700) {
                result = (1/0.85-1);
            } else{
                result = (1/0.85-1) * (1 - sin(2*M_PI*(t-700)*1/1200));
            }
            double scale = 10000; // Scale for 4 decimal places
            return std::round(result * scale) / scale;
        };
    }
    else { // no_strain
        strain_function = [](double t){
            return 0.0;
        };
    }
   
    auto runner = SimulationRunner(N_PARTICLES, SEED, strain_function);

    boost::filesystem::path executable_path = boost::filesystem::current_path();
    boost::filesystem::path root_dir = executable_path.parent_path().parent_path();

    std::cout << "Root directory: " << root_dir << std::endl;
    std::cout << "Executable path: " << executable_path << std::endl;

    std::string sequence_path;
    std::string geometry_path;
    // Use command-line arguments if provided, otherwise use default paths
    if (sequence_file.empty()) {
        sequence_path = root_dir.string() + "/sequence.yaml";
    }
    else{
        sequence_path = root_dir.string() + "/" + sequence_file;
    }

    if (geometry_file.empty()) {
        geometry_path = root_dir.string() + "/geometry_1.mat";
    }
    else{
        geometry_path = root_dir.string() + "/" + geometry_file;
    }

    //Setup sequence
    std::cout << "Sequence path: " << sequence_path << std::endl;
    std::cout << "Geometry path: " << geometry_path << std::endl;
    
    runner.setupSequence(sequence_path);
   
    //Geometry parameters
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

    // Prepare for running simulation
    parameters params;
    //Set to true if you want to save the particle positions in each time step
    params.isOutput = false;
    //Kappa parameter for the simulation - now set from command line
    params.kappa = kappa;
    //Diffusivity of the ECS
    params.D_ecs = D_ecs;
    //Diffusivity of the ICS
    params.D_ics = D_ics;
    //Number of cores for the simulation
    params.cores = cores;
    //Set to false if there is no strain (faster no need to store deformed geometries)
    params.isDeformed = isDeformed;
    //Strain step size in ms
    params.strain_step_size = strain_step_size;

    //Output parameters to file inside release folder
    std::ofstream params_file("params_job_id_" + job_id + ".txt");
    params_file << "kappa: " << kappa << std::endl;
    params_file << "seed: " << SEED << std::endl;
    params_file << "N_particles: " << N_PARTICLES << std::endl;
    params_file << "strain_type: " << strain_type << std::endl;
    params_file << "angle: " << angle << std::endl;
    params_file << "shift_block: " << shift_block << std::endl;
    params_file << "D_ecs: " << D_ecs << std::endl;
    params_file << "D_ics: " << D_ics << std::endl;
    params_file << "cores: " << cores << std::endl;
    params_file << "isDeformed: " << isDeformed << std::endl;
    params_file << "strain_step_size: " << strain_step_size << std::endl;
    params_file << "job_id: " << job_id << std::endl;
    params_file << "sequence_file: " << sequence_path << std::endl;
    params_file << "geometry_file: " << geometry_path << std::endl;
    params_file.close();
    
    // Create output filename with job ID, kappa, and seed
    std::ostringstream output_filename;
    output_filename << "sim_job" << job_id << "_" << strain_type << "_kappa" << kappa << "_seed" << SEED;
    
    //Run simulation
    runner.runSimulation(params, output_filename.str());
    
    return 0;
}
