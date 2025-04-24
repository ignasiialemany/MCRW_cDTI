//
// Created by Alemany Juvanteny, Ignasi on 23/02/2023.
//

#include "simulation.h"
#include <typeinfo>

void simulation::seedingInCuboid(const substrate &substrate, const Eigen::VectorXd &box)
{
    seeding(box, substrate);
}

void simulation::seedingInMyocytes(const substrate& substrate){
    Eigen::VectorXd box = substrate.getVoxel();
    Eigen::MatrixXd initial_positions(_particles->get_number_of_particles(), 3);
    std::mt19937 seeding(_particles->get_global_seed());
    std::uniform_real_distribution<> x_interval(box(0), box(1));
    std::uniform_real_distribution<> y_interval(box(2), box(3));
    std::uniform_real_distribution<> z_interval(box(4), box(5));

    double scale = 10000.0; // Scale for 4 decimal places

    for (int i = 0; i < _particles->get_number_of_particles(); i++)
    {
        Particle &particle = _particles->get_particle(i);
        while (true) {
            double x = x_interval(seeding);
            double y = y_interval(seeding);
            double z = z_interval(seeding);

            // Rounding the positions to 4 decimal places
            particle.position(0) = std::round(x * scale) / scale;
            particle.position(1) = std::round(y * scale) / scale;
            particle.position(2) = std::round(z * scale) / scale;

            particle.myocyte_index = substrate.searchPolygon(particle.position, 0, "global");
            if (particle.myocyte_index != -1) {
                break;
            }
        }
}
}

void simulation::seedingECS(const substrate& substrate)
{
    Eigen::VectorXd box = substrate.getVoxel();
    Eigen::MatrixXd initial_positions(_particles->get_number_of_particles(), 3);
    std::mt19937 seeding(_particles->get_global_seed());
    std::uniform_real_distribution<> x_interval(box(0), box(1));
    std::uniform_real_distribution<> y_interval(box(2), box(3));
    std::uniform_real_distribution<> z_interval(box(4), box(5));

    double scale = 10000.0; // Scale for 4 decimal places

    for (int i = 0; i < _particles->get_number_of_particles(); i++)
    {
        Particle &particle = _particles->get_particle(i);
        particle.myocyte_index = 0;
        while (particle.myocyte_index != -1) {
            double x = x_interval(seeding);
            double y = y_interval(seeding);
            double z = z_interval(seeding);

            // Rounding the positions to 4 decimal places
            particle.position(0) = std::round(x * scale) / scale;
            particle.position(1) = std::round(y * scale) / scale;
            particle.position(2) = std::round(z * scale) / scale;

            particle.myocyte_index = substrate.searchPolygon(particle.position, 0, "global");
        }
    }
}

bool simulation::seedParticlesInBox(const substrate &substrate)
{
    Eigen::VectorXd voxel_boundingbox = substrate.getVoxel();
    // Number of boxes
    try
    {
        checkBoundingBox(voxel_boundingbox);
    }
    catch (const std::runtime_error &e)
    {
        return false;
    }

    // TODO: Implement 2D logic in seeding max_point==min_point the code works for 2D,1D boxes
    // Seeding particles in bounding box
    seeding(voxel_boundingbox, substrate);
    return true;
}

void simulation::precomputeSubstrate(substrate &substrate, Eigen::VectorXd sequence_dt)
{
    if (!params.isDeformed)
    {
        Eigen::VectorXd seq_dt(2);
        seq_dt(0) = 0;
        seq_dt(1) = 0;
        substrate.preComputeSubstrate(seq_dt);
    }
    else{
        substrate.preComputeSubstrate(sequence_dt);
    }
}

boost::variant<bool, std::runtime_error> simulation::checkBoundingBox(const Eigen::VectorXd &box)
{
    int dimensions = (int)box.size() / 2;
    for (int j = 0; j < dimensions; j++)
    {
        bool isIncosistent = box(2 * j) > box(2 * j + 1);
        if (isIncosistent)
        {
            throw std::runtime_error("Bounding Box is inconsistent, min > max.");
        }
    }
    return true;
}

void simulation::seeding(const Eigen::VectorXd &box,const substrate& substrate)
{
    Eigen::MatrixXd initial_positions(_particles->get_number_of_particles(), 3);
    // Specifically for the seeding we will use the general seed as well
    std::mt19937 seeding(_particles->get_global_seed());
    std::uniform_real_distribution<> x_interval(box(0), box(1));
    std::uniform_real_distribution<> y_interval(box(2), box(3));
    std::uniform_real_distribution<> z_interval(box(4), box(5));
    for (int i = 0; i < _particles->get_number_of_particles(); i++)
    {
        Particle &particle = _particles->get_particle(i);
        particle.position(0) = x_interval(seeding);
        particle.position(1) = y_interval(seeding);
        particle.position(2) = z_interval(seeding);
        particle.myocyte_index = substrate.searchPolygon(particle.position, 0, "global");
    }
}

void simulation::writeToFile(Particle &particle)
{
    // Add particle.position, particle.myocyte_index, particle.seed, particle.index
    std::string text = std::to_string(particle.position(0)) + "," + std::to_string(particle.position(1)) + "," + std::to_string(particle.position(2)) + "," + std::to_string(particle.phase(0)) + "," + std::to_string(particle.phase(1)) + "," + std::to_string(particle.phase(2)) + "," + std::to_string(particle.myocyte_index) + "," + std::to_string(particle.flag) + "," + std::to_string(particle.exchange_time) + "\n";
    particle.buffer += text;
    if (particle.buffer.size() >= 5000) {
        particle.file.write(particle.buffer);
        particle.buffer.clear();
    }
}

void simulation::performScan(substrate &substrate_input, const sequence &sequence)
{
    if (params.isOutput)
    {
        _particles->initializeFiles();
    }

    int completed_particles = 0;
    
#pragma omp parallel num_threads(params.cores) default(none) shared(_particles, params, sequence, substrate_input,std::cout, completed_particles)
    {
        #pragma omp master // This block will be executed by only one thread (the master)
        {
            std::cout << "OpenMP is enabled. Number of threads: " << omp_get_num_threads() << std::endl;
            std::cout << "OpenMP version: " << _OPENMP << std::endl;
        }
        
        #pragma omp for schedule(static)
        for (int index_particle = 0; index_particle < _particles->get_number_of_particles(); index_particle++)
        {
            Particle &particle = _particles->get_particle(index_particle);
            particle.index = index_particle;
            /*if (particle.index <= 2020 ){
                continue;
            }*/
            try
            {
                one_walker(particle, substrate_input, sequence);
                if (particle.buffer.size() > 0) {
                    particle.file.write(particle.buffer);
                    particle.buffer.clear();
                }
            }
            catch (const std::exception &ex)
            {
#pragma omp critical
                std::cout << "Fatal error : " << index_particle << ": " << ex.what() << std::endl;
            }
#pragma omp critical
            completed_particles++;
            if (index_particle % 300 == 0)
            {
                double percentage = 100 * (double)completed_particles / (double)_particles->get_number_of_particles();
                std::cout << "Completed " << percentage << "% of particles" << std::endl;
            }
        }
    }
#pragma omp barrier
    std::cout << "Simulation finished" << std::endl;
}

void simulation::one_walker(Particle &particle, const substrate &substrate, const sequence &sequence)
{
    particle.myocyte_index = substrate.searchPolygon(particle.position, 0, "global");
    std::mt19937 rng_engine(particle.seed);
    double total_time = 0;

if (params.isOutput)
                {
                    writeToFile(particle);
                }

    
    for (int i = 0; i < sequence.dt.size(); i++)
    {
        // Get dT and dG values
        double dt_magnitude = sequence.dt(i);
        double dG_magnitude = sequence.gG(i);
        
        particle.phase = particle.phase + (dG_magnitude * dt_magnitude) * particle.position;
        // Initialize counter and flags
        bool step_success = false;
        int counter = 0;

        while (!step_success)
        {
            counter++;
            if (counter > 10)
            {
                particle.flag = 2;
                // std::cout << "Particle index COUNTER REACHED" << particle.index << " " << std::endl;
                throw std::runtime_error("Runtime error: Stepping has stopped. Particle flagged");
            }

            try
            {

                Eigen::VectorXd strain_time_array = substrate.strain_array_time;

                //Compute vectroxd with absolute values of the difference between total_time and strain_time_array
                Eigen::VectorXd strain_time_array_abs = (strain_time_array.array() - total_time).abs();

                //Find the index of the minimum value in strain_time_array_abs
                Eigen::VectorXd::Index index_min;
                strain_time_array_abs.minCoeff(&index_min);
                int index = static_cast<int>(index_min);

                //We need to find the index of the closest value to total_time in the strain_time_array
                if (!params.isDeformed)
                {
                    index = 0;
                }

                one_dt(particle, substrate, rng_engine, dt_magnitude, index, total_time);
                if (params.isOutput)
                {
                    writeToFile(particle);
                }
                step_success = true;
                total_time = total_time + dt_magnitude;
            }
            catch (const std::exception &ex)
            {
                if (typeid(ex) == typeid(std::logic_error))
                {
                    // TODO: Implement more logic errors in functions
                    particle.flag = 2;
                    std::cout << "Particle Index" << particle.index << " " << ex.what() << std::endl;
                    std::cout << "Time step :" << i << std::endl;
                    return;
                }
                else
                {
                    std::cout << "Particle Index" << particle.index << " " << ex.what() << std::endl;
                    std::string errorMessage = ex.what();
                    if (errorMessage.find("Transform") != std::string::npos)
                    {
                        // The particle local position lays very close to the block (epsilon value), just move a bit
                        particle.position = particle.position + Eigen::Vector3d::Constant(1e-6);
                        // Runtime error is in intersection, just repeat one_dt with another step
                        continue;
                    }
                    //particle.position = particle.position + Eigen::Vector3d::Constant(1e-6);
                    continue;
                }
            }
        }

    }
}

template <typename URNG>
Eigen::Vector3d simulation::getStep(URNG &rng_engine, int dimension, std::string step_type)
{
    Eigen::Vector3d step = Eigen::Vector3d::Zero(3);
    std::unordered_map<std::string, int> string_map = {
        {"constant", 0},
        {"normal", 1}};
    bool isValid = false;
    double max_step = 5;
    std::uniform_int_distribution<int> uniform_dist(0, 1);
    std::normal_distribution<double> normal(0, 1);
    while (!isValid)
    {
        switch (string_map[step_type])
        {
        case 0:
            for (int i = 0; i < dimension; i++)
            { // Maps -1 and 1 from 0 , 1
                step(i) = 2 * uniform_dist(rng_engine) - 1;
            }
            break;
        case 1:
            for (int i = 0; i < dimension; i++)
            {
                step(i) = normal(rng_engine);
            }
            break;
        }
        isValid = (step.array() * step.array()).sum() <= max_step * max_step;
    }
    return step;
}

template <typename URNG>
void simulation::one_dt(Particle &particle, const substrate &substrate, URNG &rng_engine, double dt, int index_sequence, double total_time)
{
    // Get step
    Eigen::Vector3d step = simulation::getStep(rng_engine, params.dimension, params.step_type);
    
    double dt_magnitude = dt;
    double D_coeff_old, D_coeff_new;

    // Get D coefficient
    D_coeff_old = (particle.myocyte_index != -1) ? params.D_ics : params.D_ecs;
    D_coeff_new = D_coeff_old;
    step = step * std::sqrt(2 * dt_magnitude * D_coeff_old);

    // Get variables for Transit Model (TODO: Maybe transit model as a params enumerator or class?)
    double probability_of_transit, term, D_low, l_low, p_fieremans, p_maruyama;

    int counter = 0;
    // TODO: Unify convergence epsilon in class, maybe simulation parameter?
    double convergence_eps = 1e-12;

    if (particle.strain_index_sequence != index_sequence){

        //Eigen::Vector3d reference_position = substrate.global2reference(particle.position, substrate.strain_array_value(particle.strain_index_sequence));

        double strain_value = substrate.strain_array_value(particle.strain_index_sequence);
        double angle = substrate.getAngle(particle.position, strain_value);

        double strain_value_dt = substrate.strain_array_value(index_sequence);

        Eigen::Vector3d block_centroid = substrate.getVoxelCentroid();
        Eigen::Vector3d relative_position = particle.position - block_centroid;

        Eigen::Vector3d rotated_relative_position = utility_substrate::rotate_y(relative_position, -angle);

        double x0 = rotated_relative_position(0) * std::sqrt((1 + strain_value));
        double y0 = rotated_relative_position(1) * std::sqrt((1 + strain_value));
        double z0 = rotated_relative_position(2) / (1 + strain_value);

        Eigen::Vector3d next_relative_position;

        next_relative_position(0) = x0 / std::sqrt((1 + strain_value_dt));
        next_relative_position(1) = y0 / std::sqrt((1 + strain_value_dt));
        next_relative_position(2) = (1 + strain_value_dt) * z0;

        Eigen::Vector3d displacement = next_relative_position - rotated_relative_position;

        Eigen::Vector3d new_disp = utility_substrate::rotate_y(displacement, angle);

        particle.position = particle.position + new_disp;
        particle.strain_index_sequence = index_sequence;
    }
    
    double strain_z = substrate.strain_array_value(index_sequence);
    //transform_info transform_data = substrate.getLocalFromGlobal(particle.position,strain_z);
    //Eigen::Vector3d local_position = transform_data.local_position;
    //Eigen::Vector3d local_step = utility_substrate::rotate_y(step, -transform_data.angle);

    //We iterate until the step is less than the convergence
    while (step.norm() > convergence_eps)
    {
        //We will perform global to local, deal with intersections and then local to global againa

        //GLOBAL TO LOCAL from particle.position and step
        transform_info transform_data = substrate.getLocalFromGlobal(particle.position,strain_z);
        Eigen::Vector3d local_position = transform_data.local_position;
        Eigen::Vector3d local_step = utility_substrate::rotate_y(step, -transform_data.angle);

        //INIT REMAINING STEP 
        Eigen::Vector3d local_normalized_step = local_step / local_step.norm();
        Eigen::Vector3d remaining_step;

        D_coeff_old = D_coeff_new;
        auto index_polygon = substrate.searchPolygon(local_position, index_sequence);
        //Update D_coeff_old as expected
        if (particle.myocyte_index!=-1){
            if (index_polygon==-1){
                throw std::logic_error("FLAGGED - ICS to ECS");
               D_coeff_old = params.D_ecs;
            }
            D_coeff_old = params.D_ics;
        }
        else{
            if (index_polygon!=-1){
                throw std::logic_error("FLAGGED - ECS to ICS");
                D_coeff_old = params.D_ics;
            }
            D_coeff_old = params.D_ecs;
        }
        

        counter++;
        //  Throw logic_error if we we try to intersect more than 50 times
        if (counter > 50)
        {
            throw std::logic_error("FLAGGED - Stepping error, stopped while loop in one_dt after trying 50 times");
        }
        
        // Get intersection data
        boost::optional<std::tuple<int, double, Eigen::Vector3d>> intersection_data = substrate.intersectPolygon(local_position, local_step, index_sequence);

        // If intersection is type bool and false
        if (intersection_data)
        {
            // Get intersection data
            int polygon_index = std::get<0>(*intersection_data);
            double distance_to_intersection = std::get<1>(*intersection_data);
            Eigen::Vector3d normal = std::get<2>(*intersection_data);

            // Obtain the remaining step and the step to the intersection
            double remaining_value = local_step.norm() - distance_to_intersection;
            remaining_step = local_normalized_step * remaining_value;
            Eigen::Vector3d step_to_intersection = local_normalized_step * distance_to_intersection;

            // dot product between normal and normalized step
            double dot_product = normal.dot(local_normalized_step);
            double cos_angle = dot_product > 0 ? dot_product : -dot_product;
            double distance_to_intersection_projected = distance_to_intersection * cos_angle;

            // TODO: Implement transit models, for now we consider the hybrid model

            if (D_coeff_old == params.D_ecs)
            {
                D_low = params.D_ics;
                l_low = distance_to_intersection_projected * std::sqrt(D_low / D_coeff_old);
                term = 2 * l_low * params.kappa / D_low;
                p_fieremans = term / (term + 1);
                p_maruyama = std::sqrt(params.D_ics / D_coeff_old) > 1 ? 1 : std::sqrt(params.D_ics / D_coeff_old); // equivalent to matlab min(1,sqrt(substrate.D_i/D_old))
                probability_of_transit = p_fieremans * p_maruyama;
            }
            else
            {
                term = 2 * distance_to_intersection_projected * params.kappa / params.D_ics;
                p_fieremans = term / (term + 1);
                // p_maruyama = min(1,sqrt(D_new/D_old)); (in matlab)
                probability_of_transit = p_fieremans;
            }

            std::uniform_real_distribution<> transit(0, 1);
            double U;
            if (typeid(rng_engine) == typeid(oneGenerator))
            {
                U = 0.5; // Hardcoded value for U when rng_engine is of type oneGenerator
            }
            else
            {
                U = transit(rng_engine);
            }

            // Particle is crossing
            if (U < probability_of_transit)
            {
                if(params.isOutput){
                    particle.exchange_time = total_time - particle.exchange_time;
                    writeToFile(particle);
                    particle.exchange_time = total_time;
                }
                // If the particle is in ECS update to myocyte index otherwise update to -1
                if (particle.myocyte_index == -1)
                {
                    particle.myocyte_index = polygon_index;
                    D_coeff_new = params.D_ics;
                }
                else
                {
                    particle.myocyte_index = -1;
                    D_coeff_new = params.D_ecs;
                }
                // Update remaining step to account for the change in D
                remaining_step = remaining_step * std::sqrt(D_coeff_new / D_coeff_old);
            }
            else
            {

                // Compute the reflection based on the remaining normalized step
                Eigen::Vector3d remaining_step_normalized = remaining_step / remaining_step.norm();
                Eigen::Vector3d reflected_step = remaining_step_normalized - 2 * normal * normal.dot(remaining_step_normalized);
                // Reflect the remaining_step
                remaining_step = reflected_step * remaining_step.norm();

                if (substrate.searchPolygon(local_position, index_sequence)!=-1){
                int poly = substrate.searchPolygon(local_position, index_sequence);
                substrate.searchPolygon(local_position, index_sequence);
                D_coeff_old = params.D_ics;
            }

            }

            local_position = local_position + step_to_intersection;
            local_position = local_position + remaining_step * 1e-6;
            remaining_step = remaining_step * (1 - 1e-6);
        }
        else
        {
            //TODO: Implement boundary conditions, "reflective" and "periodic"(currently is the one that is implemented if it's not none)
            if(substrate.getBoundaryType() == "none")
            {
                if (substrate.searchPolygon(local_position, index_sequence)!=-1){
                int poly = substrate.searchPolygon(local_position, index_sequence);
                substrate.searchPolygon(local_position, index_sequence);
                D_coeff_old = params.D_ics;
                }

                local_position = local_position + local_step;
                remaining_step = Eigen::Vector3d::Zero(3);
                local_step = remaining_step;
                continue;
            }

            // Get intersection with block
            boost::optional<std::tuple<int, double, Eigen::Vector3d>> intersection_block = substrate.intersectionBlock(local_position, local_step, index_sequence);

            // If it intersects with block
            if (intersection_block)
            {

                //std::cout << "CROSSED BLOCK" << std::endl;

                int block_face = std::get<0>(*intersection_block);
                double distance = std::get<1>(*intersection_block);
                Eigen::Vector3d normal = std::get<2>(*intersection_block);

                //To make sure it lands on the other block
                // Compute remaining_step and step_to_intersection
                remaining_step = local_normalized_step * (local_step.norm() - distance);
                Eigen::Vector3d step_to_intersection = local_normalized_step * distance;
                
                step_to_intersection = step_to_intersection  + local_normalized_step*1e-6;
                remaining_step = remaining_step -  local_normalized_step*1e-6;

                if (substrate.getBoundaryType() == "reflective")
                {
                     Eigen::Vector3d reflected_step_normalized = local_normalized_step - 2 * normal * normal.dot(local_normalized_step);
                     remaining_step = reflected_step_normalized * remaining_step.norm();
                     local_position = local_position + reflected_step_normalized * 1e-6;
                    remaining_step = remaining_step * (1 - 1e-6);
                }
               
                local_position = local_position + step_to_intersection;
                int possible_myo = substrate.searchPolygon(local_position, index_sequence);

                if (possible_myo != -1)
                {
                    throw std::runtime_error("Particle has crossed block and directly to cardiomyocyte");
                }
            }
            else
            {
                local_position = local_position + local_step;
                remaining_step = Eigen::Vector3d::Zero(3);
            }
        }
        local_step = remaining_step;
        //Transform local_step and local_position to step and particle.position back again
        step = utility_substrate::rotate_y(local_step, transform_data.angle);
        particle.position = substrate.getGlobalFromLocal(local_position, transform_data.iX, transform_data.iY, transform_data.iZ, strain_z);
    }
}

template void simulation::one_dt<oneGenerator>(Particle &particle, const substrate &substrate, oneGenerator &rng_engine, double dt, int index_sequence, double total_time);
template void simulation::one_dt<std::mt19937>(Particle &particle, const substrate &substrate, std::mt19937 &rng_engine, double dt, int index_sequence, double total_time);
template Eigen::Vector3d simulation::getStep<oneGenerator>(oneGenerator &rng_engine, int dimension, std::string step_type);
template Eigen::Vector3d simulation::getStep<std::mt19937>(std::mt19937 &rng_engine, int dimension, std::string step_type);
