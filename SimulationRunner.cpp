// SimulationRunner.cpp
#include "SimulationRunner.h"


SimulationRunner::SimulationRunner(int N_particles, int seed, std::function<double(double)> strain_function) 
    : N_particles_(N_particles),
      seed_(seed),
      strain_function_(strain_function)
{
    simulation_ptr_ = std::make_unique<simulation>(N_particles_, seed_);    
}

void SimulationRunner::setupGeometryMatfile(double angle, bool shift_block, std::string geometry_path, Eigen::VectorXd voxel, Eigen::VectorXd block_size){
    //Setup transform - we set to False identity, identity is used when there is no transform from global to local and vice versa
    transform_ptr_ = std::make_unique<transform>(angle, shift_block, false);
    transform_ptr_->set_block(block_size[0], block_size[1], block_size[2]);

    //Setup substrate using geometry_path (.mat file) from previous papers
    auto data = utility_substrate::read_mat_file(geometry_path);
    substrate_ptr_ = std::make_unique<substrate>(data.first, data.second, *transform_ptr_, strain_function_);
    substrate_ptr_->setVoxel(voxel);
}

void SimulationRunner::setupGeometryOffFiles(double angle, bool shift_block, std::vector<std::string> geometry_paths, Eigen::VectorXd voxel, Eigen::VectorXd block_size){
    
    //Setup transform - we set to False identity, identity is used when there is no transform from global to local and vice versa
    transform_ptr_ = std::make_unique<transform>(angle, shift_block, false);
    transform_ptr_->set_block(block_size[0], block_size[1], block_size[2]);

    //Setup substrate using geometry_path (.off files)
    substrate_ptr_ = std::make_unique<substrate>(geometry_paths, *transform_ptr_, strain_function_);
    substrate_ptr_->setVoxel(voxel);
}

void SimulationRunner::setupSequence(std::string seq_path){
    sequence_ptr_ = std::make_unique<sequence>(seq_path);
    sequence_ptr_->create();
    //Compute b-value and print it
    Eigen::VectorXd cumulative_dt = computeCumulativeArray(sequence_ptr_->dt);
    auto cumulative_integral = cumulativeIntegral(cumulative_dt, sequence_ptr_->gG);
    Eigen::VectorXd squared_cumulative_integral = cumulative_integral.array().square();
    auto bvalue = computeIntegral(cumulative_dt, squared_cumulative_integral);
    std::cout << "--------------------------------" << std::endl;
    std::cout << "SEQUENCE BVALUE: " << bvalue << std::endl;
    std::cout << "SEQUENCE TYPE: " << sequence_ptr_->parameters.type << std::endl;
    std::cout << "N TIMESTEPS: " << sequence_ptr_->dt.size() << std::endl;
    std::cout << "N PARTICLES: " << N_particles_ << std::endl;
    std::cout << "--------------------------------" << std::endl;
}

Eigen::VectorXd SimulationRunner::createStrainTimeArray(double step_size){
    double total_time = sequence_ptr_->dt.sum();
    int numElements = static_cast<int>(std::floor(total_time / step_size)) + 1;
    auto strain_array_time = Eigen::VectorXd(numElements);
    for (int i = 0; i < numElements - 1; ++i) {
        strain_array_time(i) = step_size;
    }
    strain_array_time(numElements - 1) = total_time - strain_array_time.head(numElements - 1).sum();
    return strain_array_time;
}

void SimulationRunner::runSimulation(parameters params, std::string file_name){

    //Set parameters
    simulation_ptr_->params = params;

    //Compute strain time array for the simulation
    auto strain_array_time = createStrainTimeArray(params.strain_step_size);
    auto start = std::chrono::high_resolution_clock::now();
    //Precompute substrate
    simulation_ptr_->precomputeSubstrate(*substrate_ptr_, strain_array_time);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Precomputed substrate took: " << elapsed.count() << " s\n";
    std::cout << "--------------------------------" << std::endl;
    //Seed particles
    simulation_ptr_->seedParticlesInBox(*substrate_ptr_);
    std::cout << "--------------------------------" << std::endl;
    //Write CSV file with initial values
    simulation_ptr_->writeParticlesState(file_name + "_init.csv", *substrate_ptr_);
    std::cout << "--------------------------------" << std::endl;
    //Run simulation
    auto start_sim = std::chrono::high_resolution_clock::now();
    simulation_ptr_->performScan(*substrate_ptr_, *sequence_ptr_);
    auto end_sim = std::chrono::high_resolution_clock::now();   
    std::chrono::duration<double> elapsed_sim = end_sim - start_sim;
    std::cout << "Simulation took: " << elapsed_sim.count() << " s\n";
    std::cout << "--------------------------------" << std::endl;
    //Write CSV file with final values
    simulation_ptr_->writeParticlesState(file_name + "_final.csv", *substrate_ptr_);
    std::cout << "--------------------------------" << std::endl;
}


// Static utility functions
Eigen::VectorXd SimulationRunner::computeCumulativeArray(const Eigen::VectorXd& x) {
    Eigen::VectorXd cumulative(x.size());
    cumulative(0) = x(0);
    for (int i = 1; i < x.size(); ++i) {
        cumulative(i) = cumulative(i - 1) + x(i);
    }
    return cumulative;
}

Eigen::VectorXd SimulationRunner::cumulativeIntegral(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
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

double SimulationRunner::computeIntegral(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
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