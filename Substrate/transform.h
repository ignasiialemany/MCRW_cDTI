//
// Created by Alemany Juvanteny, Ignasi on 26/02/2023.
//

#ifndef INC_3DRANDOMWALK_TRANSFORM_H
#define INC_3DRANDOMWALK_TRANSFORM_H

#include <Eigen/Dense>
#include "utility_substrate.h"
#include <cmath>
#include <CGAL/Bbox_3.h>
#include "../Geometry/polyhedronSet.h"

struct transform_info{
    Eigen::Vector3d local_position;
    //TODO: Might delete angle 
    double angle;
    int iX,iY,iZ;
};


class transform {

public:
    transform()=default;
    transform(double rot_in_y, bool shift_block, bool isIdentity): deg_rot_per_m_in_Y(rot_in_y), shift_block(shift_block), isIdentity(isIdentity){};

    void precomputeStrainedBlock();
    transform_info global2local(const Eigen::Vector3d &global_position, double strain_z, const Eigen::Vector3d& centroid_voxel) const;
    Eigen::Vector3d local2global(const Eigen::Vector3d &local_position, int iX, int iY, int iZ,  double strain_z, const Eigen::Vector3d& centroid_voxel) const;
    void set_block(double dx, double dy, double dz, double minY = -15000. , double maxY = 15000.,  double minX=-15000., double maxX=15000., double minZ=-15000., double maxZ=15000.);
    Eigen::Vector3d get_block_size()const {return Eigen::Vector3d(dx,dy,dz);};
    void precomputeTransform(Eigen::VectorXd &sequence_dt, std::function<double(double)> &strain);
    bool isTransformIdentity() const {return isIdentity;};
    Eigen::Vector3d get_block_centroid() const {return _centroid_block;};
    double getAngle(const Eigen::Vector3d &global_position, double strain_z, const Eigen::Vector3d &centroid_voxel) const;
    const polygon& getBlockCurrentTime(int sequence_index) const {return _block.getPolygon(sequence_index);};
    //Eigen::Vector3d global2reference(const Eigen::Vector3d &point, double strain_z, const Eigen::Vector3d& centroid_voxel) const;
    //Eigen::Vector3d reference2global(const Eigen::Vector3d &point, int iZ, double angle, double strain_z) const;
 Eigen::MatrixXd y_slice_minmax;
    Eigen::MatrixXd x_slice_minmax;
    Eigen::MatrixXd z_slice_minmax;

private:

    void create_block(Kernel::Point_3 min_point, Kernel::Point_3 max_point);
    double deg_rot_per_m_in_Y = 0.01;
    void calculate_minmax_slices(Eigen::MatrixXd& slice_minmax, double minrange, double maxrange);
   
    bool shift_block = true;
    double dx,dy,dz;
    bool isIdentity;
    std::vector<Kernel::Triangle_3> triangles;
    polyhedronSet _block;
    Eigen::Vector3d _centroid_block;
};


#endif //INC_3DRANDOMWALK_TRANSFORM_H
