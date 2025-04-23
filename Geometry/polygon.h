//
// Created by Alemany Juvanteny, Ignasi on 27/02/2023.
//

#ifndef INC_3DRANDOMWALK_POLYGON_H
#define INC_3DRANDOMWALK_POLYGON_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Ray_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/intersection_3.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Point_3.h>
#include <Eigen/Dense>
#include <cassert>
#include <vector>
#include <boost/optional.hpp>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk_cost.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/property_map.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/IO/PLY.h>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include <CGAL/Object.h>
#include <memory>
#include <fstream>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Surface_mesh/IO/PLY.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Side_of_triangle_mesh.h>
//typedef CGAL::Simple_cartesian<long double> Kernel;
typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::HalfedgeDS HalfedgeDS;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;

// Define the AABB traits
typedef std::vector<Kernel::Triangle_3>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<Kernel, Iterator> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree_AABB;
typedef CGAL::Exact_predicates_exact_constructions_kernel ExactKernel;

class polygon
{
public:
    polygon() = default;
    //Add deconstructor
    polygon(const Eigen::MatrixXd &vertices_input, const Eigen::MatrixXd &faces_input);
    polygon(Polyhedron &poly);
    polygon(const polygon& other);
    polygon& operator=(const polygon& other);


    double computeVolume();
    double computeSurface();
    Polyhedron getPolyhedron() { return _poly; };
    Polyhedron getBbox() { return _bbox; };
    CGAL::Bbox_3 getSolidBbox() const { return _solid_bbox; };
    
    boost::optional<std::pair<int, double>> intersection(const Eigen::Vector3d &point, const Eigen::Vector3d &step) const;
    //boost::optional<std::pair<int, double>> intersectionSecond(const Eigen::Vector3d &point, const Eigen::Vector3d &step) const;

    //boost::optional<std::pair<int, double>> intersectionSecond(const Eigen::Vector3d &point, const Eigen::Vector3d &step) const;
    //boost::optional<std::pair<int, double>> intersectionForBlock(const Eigen::Vector3d &point, const Eigen::Vector3d &step) const;

    // bool containsPoint(const Eigen::Vector3d &point);
    bool containsPoint(const Eigen::Vector3d &point) const;
    Eigen::Vector3d getNormalVector(int index_face) const;
    std::vector<CGAL::Point_3<Kernel>> getVertices() const { return _vertices; };
    std::vector<std::vector<std::size_t>> getFaces() const { return _faces; };
    bool isPolygonClosed();
    Mesh getMesh() { return mesh; };
    void updatePolygon(double strain, Eigen::Vector3d& centroid);
    Eigen::Vector3d getCentroid() const { return centroid; };
        std::vector<Kernel::Triangle_3> triangle_faces;


private:
    // TODO: Might delete some of these variables, probably _vertices and _faces
    Polyhedron _poly;
    Polyhedron _bbox;
    CGAL::Bbox_3 _solid_bbox;
    CGAL::Side_of_triangle_mesh<Polyhedron,Kernel> _side_of_triangle_mesh;
    Mesh mesh;  
    void createBbox(Kernel::Point_3 min_point, Kernel::Point_3 max_point);
    void createPolygon(const Eigen::MatrixXd &vertices, const Eigen::MatrixXd &faces);
    std::vector<CGAL::Point_3<Kernel>> _vertices;
    std::vector<std::vector<std::size_t>> _faces;
    std::unique_ptr<Tree_AABB> _AABBtree;
    Eigen::Vector3d centroid;
};

#endif // INC_3DRANDOMWALK_POLYGON_H
