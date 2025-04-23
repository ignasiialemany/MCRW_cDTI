//
// Created by Alemany Juvanteny, Ignasi on 26/02/2023.
//

#include "transform.h"
#include <unordered_set>

double transform::getAngle(const Eigen::Vector3d &global_position, double strain_z,const Eigen::Vector3d &centroid_voxel) const
{
    double strain_x = std::sqrt(1/(1+strain_z)) - 1;
    double strain_y = std::sqrt(1/(1+strain_z)) - 1;
    double dx_strained = dx * (1 + strain_x);
    double dy_strained = dy * (1 + strain_y);
    double dz_strained = dz * (1 + strain_z);

    if (isIdentity)
    {
        return 0;
    }
    else
    {
        //calculate new slices with strained values
        Eigen::MatrixXd strained_z_slice_minmax;
        Eigen::MatrixXd strained_y_slice_minmax;
        Eigen::MatrixXd strained_x_slice_minmax;
        //Mltiply z_slice_minmax Eigen::MatrixXd by (1+strain_z)

        strained_z_slice_minmax.resize(2, z_slice_minmax.cols());
        strained_z_slice_minmax.row(0) = z_slice_minmax.row(0) * (1 + strain_z);
        strained_z_slice_minmax.row(1) = z_slice_minmax.row(1) * (1 + strain_z);

        strained_y_slice_minmax.resize(2, y_slice_minmax.cols());
        strained_y_slice_minmax.row(0) = y_slice_minmax.row(0) * (1 + strain_y);
        strained_y_slice_minmax.row(1) = y_slice_minmax.row(1) * (1 + strain_y);

        strained_x_slice_minmax.resize(2, x_slice_minmax.cols());
        strained_x_slice_minmax.row(0) = x_slice_minmax.row(0) * (1 + strain_x);
        strained_x_slice_minmax.row(1) = x_slice_minmax.row(1) * (1 + strain_x);
       
        double disp_x = centroid_voxel(0)*(1 + strain_x) - centroid_voxel(0);
        double disp_y = centroid_voxel(1)*(1 + strain_y) - centroid_voxel(1);
        double disp_z = centroid_voxel(2)*(1 + strain_z) - centroid_voxel(2);

        strained_z_slice_minmax.row(0) = strained_z_slice_minmax.row(0).array()  - disp_z;
        strained_z_slice_minmax.row(1) = strained_z_slice_minmax.row(1).array()  - disp_z;
        strained_y_slice_minmax.row(0) = strained_y_slice_minmax.row(0).array()  - disp_y;
        strained_y_slice_minmax.row(1) = strained_y_slice_minmax.row(1).array()  - disp_y;
        strained_x_slice_minmax.row(0) = strained_x_slice_minmax.row(0).array()  - disp_x;
        strained_x_slice_minmax.row(1) = strained_x_slice_minmax.row(1).array()  - disp_x;
        
        int iY = utility_substrate::find_index_slice(global_position(1), strained_y_slice_minmax);
        double y_slice = y_slice_minmax(0, iY);
        //double y_slice = utility_substrate::find_yslice(global_position(1), strained_y_slice_minmax);
        double angle = (deg_rot_per_m_in_Y * M_PI / 180.) * (y_slice); // This angle is in rad now
        double angle_rounded= std::round(angle * 100) / 100;
        return angle_rounded;
    }
}


transform_info transform::global2local(const Eigen::Vector3d &global_position, double strain_z, const Eigen::Vector3d& centroid_voxel) const
{
    transform_info output;

    double strain_x = std::sqrt(1/(1+strain_z)) - 1;
    double strain_y = std::sqrt(1/(1+strain_z)) - 1;
    double dx_strained = dx * (1 + strain_x);
    double dy_strained = dy * (1 + strain_y);
    double dz_strained = dz * (1 + strain_z);

    if (isIdentity)
    {
        output.angle = 0;
        output.local_position = global_position;
        return output;
    }
    else
    {
        //calculate new slices with strained values
        Eigen::MatrixXd strained_z_slice_minmax;
        Eigen::MatrixXd strained_y_slice_minmax;
        Eigen::MatrixXd strained_x_slice_minmax;
        //Mltiply z_slice_minmax Eigen::MatrixXd by (1+strain_z)

        strained_z_slice_minmax.resize(2, z_slice_minmax.cols());
        strained_z_slice_minmax.row(0) = z_slice_minmax.row(0) * (1 + strain_z);
        strained_z_slice_minmax.row(1) = z_slice_minmax.row(1) * (1 + strain_z);

        strained_y_slice_minmax.resize(2, y_slice_minmax.cols());
        strained_y_slice_minmax.row(0) = y_slice_minmax.row(0) * (1 + strain_y);
        strained_y_slice_minmax.row(1) = y_slice_minmax.row(1) * (1 + strain_y);

        strained_x_slice_minmax.resize(2, x_slice_minmax.cols());
        strained_x_slice_minmax.row(0) = x_slice_minmax.row(0) * (1 + strain_x);
        strained_x_slice_minmax.row(1) = x_slice_minmax.row(1) * (1 + strain_x);
       
        double disp_x = centroid_voxel(0)*(1 + strain_x) - centroid_voxel(0);
        double disp_y = centroid_voxel(1)*(1 + strain_y) - centroid_voxel(1);
        double disp_z = centroid_voxel(2)*(1 + strain_z) - centroid_voxel(2);

        strained_z_slice_minmax.row(0) = strained_z_slice_minmax.row(0).array()  - disp_z;
        strained_z_slice_minmax.row(1) = strained_z_slice_minmax.row(1).array()  - disp_z;
        strained_y_slice_minmax.row(0) = strained_y_slice_minmax.row(0).array()  - disp_y;
        strained_y_slice_minmax.row(1) = strained_y_slice_minmax.row(1).array()  - disp_y;
        strained_x_slice_minmax.row(0) = strained_x_slice_minmax.row(0).array()  - disp_x;
        strained_x_slice_minmax.row(1) = strained_x_slice_minmax.row(1).array()  - disp_x;

        output.iY = utility_substrate::find_index_slice(global_position(1), strained_y_slice_minmax);

        double y_slice = y_slice_minmax(0, output.iY);
        //double y_slice = utility_substrate::find_yslice(global_position(1), strained_y_slice_minmax);
        output.angle = (deg_rot_per_m_in_Y * M_PI / 180.) * (y_slice); // This angle is in rad now
        output.angle = std::round(output.angle * 100) / 100;
        
        Eigen::Vector3d relative_pos = global_position - centroid_voxel;
        
        Eigen::Vector3d rel_position_rotated = utility_substrate::rotate_y(relative_pos, -output.angle);

        Eigen::Vector3d position_rotated = rel_position_rotated + centroid_voxel;

        //find the slice and index of the position
        output.iZ = utility_substrate::find_index_slice(position_rotated(2), strained_z_slice_minmax);
        
        if (shift_block){
        if (output.iZ%2 == 1){
            position_rotated(0) = position_rotated(0) - dx_strained/2;
        }
        }
        
        output.iX = utility_substrate::find_index_slice(position_rotated(0), strained_x_slice_minmax);
        
        //Let's compute the centroid of the block they are in from the strained_x_slice_minmax.col(output.iX)
        double centroid_x  = (strained_x_slice_minmax(0, output.iX) + strained_x_slice_minmax(1, output.iX))/2;
        double centroid_y  = (strained_y_slice_minmax(0, output.iY) + strained_y_slice_minmax(1, output.iY))/2;
        double centroid_z  = (strained_z_slice_minmax(0, output.iZ) + strained_z_slice_minmax(1, output.iZ))/2;
        Eigen::Vector3d centroid_block = Eigen::Vector3d(centroid_x, centroid_y, centroid_z);

        // Translation offset
        output.local_position = position_rotated - centroid_block;
        
        // throw error if local_position is not within block but consider epsilon
        if (std::abs(output.local_position(0))-dx_strained/2 > 1e-5 || std::abs(output.local_position(1))-dy_strained/2 > 1e-5 || std::abs(output.local_position(2))-dz_strained/2 > 1e-5)
        {
            //print out the eps diff between the local_position and the block edges
            std::cout << "eps_x" << std::abs(output.local_position(0)) - dx_strained/2 << std::endl;
            std::cout << "eps_y" << std::abs(output.local_position(1)) - dy_strained/2 << std::endl;
            std::cout << "eps_z" << std::abs(output.local_position(2)) - dz_strained/2 << std::endl;
            

            std::cout << "Local position: " << output.local_position << std::endl;
            throw std::runtime_error("Transform::global2local -> Local position is not within block");
        }
       
         
        return output;
    }
}
Eigen::Vector3d transform::local2global(const Eigen::Vector3d &local_position, int iX, int iY, int iZ, double strain_z, const Eigen::Vector3d& centroid_voxel) const
{
    if (isIdentity)
    {
        return local_position;
    }

    double strain_x = std::sqrt(1/(1+strain_z)) - 1;
    double strain_y = std::sqrt(1/(1+strain_z)) - 1;

    double dz_strained = dz * (1 + strain_z);
    double dy_strained = dy * (1 + strain_y);
    double dx_strained = dx * (1 + strain_x);

    Eigen::MatrixXd strained_z_slice_minmax;
    Eigen::MatrixXd strained_y_slice_minmax;
    Eigen::MatrixXd strained_x_slice_minmax;
    
    strained_z_slice_minmax.resize(2, z_slice_minmax.cols());
    strained_z_slice_minmax.row(0) = z_slice_minmax.row(0) * (1 + strain_z);
    strained_z_slice_minmax.row(1) = z_slice_minmax.row(1) * (1 + strain_z);

    strained_y_slice_minmax.resize(2, y_slice_minmax.cols());
    strained_y_slice_minmax.row(0) = y_slice_minmax.row(0) * (1 + strain_y);
    strained_y_slice_minmax.row(1) = y_slice_minmax.row(1) * (1 + strain_y);

    strained_x_slice_minmax.resize(2, x_slice_minmax.cols());
    strained_x_slice_minmax.row(0) = x_slice_minmax.row(0) * (1 + strain_x);
    strained_x_slice_minmax.row(1) = x_slice_minmax.row(1) * (1 + strain_x);

    double disp_x = centroid_voxel(0)*(1 + strain_x) - centroid_voxel(0);
    double disp_y = centroid_voxel(1)*(1 + strain_y) - centroid_voxel(1);
    double disp_z = centroid_voxel(2)*(1 + strain_z) - centroid_voxel(2);

    strained_z_slice_minmax.row(0) = strained_z_slice_minmax.row(0).array()  - disp_z;
    strained_z_slice_minmax.row(1) = strained_z_slice_minmax.row(1).array()  - disp_z;
    strained_y_slice_minmax.row(0) = strained_y_slice_minmax.row(0).array()  - disp_y;
    strained_y_slice_minmax.row(1) = strained_y_slice_minmax.row(1).array()  - disp_y;
    strained_x_slice_minmax.row(0) = strained_x_slice_minmax.row(0).array()  - disp_x;
    strained_x_slice_minmax.row(1) = strained_x_slice_minmax.row(1).array()  - disp_x;

    double centroid_x  = (strained_x_slice_minmax(0, iX) + strained_x_slice_minmax(1, iX))/2;
    double centroid_y  = (strained_y_slice_minmax(0, iY) + strained_y_slice_minmax(1, iY))/2;
    double centroid_z  = (strained_z_slice_minmax(0, iZ) + strained_z_slice_minmax(1, iZ))/2;
    Eigen::Vector3d centroid_block = Eigen::Vector3d(centroid_x, centroid_y, centroid_z);

    Eigen::Vector3d position_global = local_position + centroid_block;

    if (shift_block){
    if (iZ%2 == 1){
            position_global(0) = position_global(0) + dx_strained/2;
    }
    }

    double y_slice = y_slice_minmax(0, iY);
    //double y_slice = utility_substrate::find_yslice(position_global(1), strained_y_slice_minmax);
    double angle = (deg_rot_per_m_in_Y * M_PI / 180.) * y_slice; // This angle is in rad now
    angle = std::round(angle * 100) / 100; // Round to 2 decimal places

            //Correct for the rotation around the origin as it should be through the centroid
    Eigen::Vector3d relative_pos = position_global - centroid_voxel;
        
    Eigen::Vector3d rel_position_rotated = utility_substrate::rotate_y(relative_pos, angle);

    Eigen::Vector3d position_rotated = rel_position_rotated + centroid_voxel;
        

    return position_rotated;
}
void transform::set_block(double dx, double dy, double dz, double minY, double maxY, double minX, double maxX, double minZ, double maxZ)
{
    this->dx = dx;
    this->dy = dy;
    this->dz = dz;
    Kernel::Point_3 min_point = Kernel::Point_3(0, 0, 0);
    Kernel::Point_3 max_point = Kernel::Point_3(dx, dy, dz);

    //Compute centroid from min_point and max_point
    _centroid_block = Eigen::Vector3d(dx/2, dy/2, dz/2);
    create_block(min_point, max_point);

    std::unordered_set<double> y_minvals_set;
    for (double y = 0.0; y <= maxY; y += dy)
    {
        y_minvals_set.insert(y);
    }
    for (double y = -dy; y >= minY; y -= dy)
    {
        y_minvals_set.insert(y);
    }

    Eigen::VectorXd y_minvals(y_minvals_set.size());
    std::copy(y_minvals_set.begin(), y_minvals_set.end(), y_minvals.data());
    std::sort(y_minvals.data(), y_minvals.data() + y_minvals.size());

    this->y_slice_minmax.resize(2, y_minvals.size() - 1);
    this->y_slice_minmax.row(0) = y_minvals.segment(0, y_minvals.size() - 1);
    this->y_slice_minmax.row(1) = y_minvals.tail(y_minvals.size() - 1);


    std::cout << "y_slice_minmax" << std::endl;
    std::cout << y_slice_minmax << std::endl;

    std::unordered_set<double> x_minvals_set;
    for (double x = 0.0; x <= maxX; x += dx)
    {
        x_minvals_set.insert(x);
    }
    for (double x = -dx; x >= minX; x -= dx)
    {
        x_minvals_set.insert(x);
    }

    Eigen::VectorXd x_minvals(x_minvals_set.size());
    std::copy(x_minvals_set.begin(), x_minvals_set.end(), x_minvals.data());
    std::sort(x_minvals.data(), x_minvals.data() + x_minvals.size());

    this->x_slice_minmax.resize(2, x_minvals.size() - 1);
    this->x_slice_minmax.row(0) = x_minvals.segment(0, x_minvals.size() - 1);
    this->x_slice_minmax.row(1) = x_minvals.tail(x_minvals.size() - 1);
    
    std::unordered_set<double> z_minvals_set;
    for (double z = 0.0; z <= maxZ; z += dz)
    {
        z_minvals_set.insert(z);
    }
    for (double z = -dz; z >= minZ; z -= dz)
    {
        z_minvals_set.insert(z);
    }

    Eigen::VectorXd z_minvals(z_minvals_set.size());
    std::copy(z_minvals_set.begin(), z_minvals_set.end(), z_minvals.data());
    std::sort(z_minvals.data(), z_minvals.data() + z_minvals.size());

    this->z_slice_minmax.resize(2, z_minvals.size() - 1);
    this->z_slice_minmax.row(0) = z_minvals.segment(0, z_minvals.size() - 1);
    this->z_slice_minmax.row(1) = z_minvals.tail(z_minvals.size() - 1);
}

void transform::create_block(Kernel::Point_3 min_point, Kernel::Point_3 max_point)
{
    Eigen::MatrixXd ver(8,3);
    Eigen::MatrixXd faces(12,3);
    ver << min_point.x(), min_point.y(), min_point.z(),
                max_point.x(), min_point.y(), min_point.z(),
                max_point.x(), max_point.y(), min_point.z(),
                min_point.x(), max_point.y(), min_point.z(),
                min_point.x(), min_point.y(), max_point.z(),
                max_point.x(), min_point.y(), max_point.z(),
                max_point.x(), max_point.y(), max_point.z(),
                min_point.x(), max_point.y(), max_point.z();
    faces << 0, 1, 2,
             0, 2, 3,
             4, 6, 5,
             4, 7, 6,
             0, 4, 5,
             0, 5, 1,
             1, 5, 6,
             1, 6, 2,
             2, 6, 7,
             2, 7, 3,
             3, 7, 4,
             3, 4, 0;
    ;
    faces = faces.array() + 1;
    _block = polyhedronSet(ver, faces);    
}


void transform::precomputeTransform(Eigen::VectorXd &sequence_dt, std::function<double(double)> &strain){
    _block.precomputePolygons(_centroid_block,sequence_dt,strain);
}
