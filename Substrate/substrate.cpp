//
// Created by Alemany Juvanteny, Ignasi on 26/02/2023.
//

#include "substrate.h"
#include <omp.h>
#include <chrono>

substrate::substrate(std::vector<Eigen::MatrixXd> &myo_vertices, std::vector<Eigen::MatrixXd> &myo_faces, transform &t, std::function<double(double)> strain_f)
{
    _transform = std::move(t);
    _strain = strain_f;
    _myocytes.reserve(myo_vertices.size());
    int i=0;
    for (auto &myo_vertices_i : myo_vertices)
    {
        polyhedronSet myo(myo_vertices_i, myo_faces[i]);
        Eigen::Vector3d centroid = myo.getUnstrainedCentroid();
        Kernel::Point_3 centroid_cgal(centroid(0), centroid(1), centroid(2));
        _points.emplace_back(centroid_cgal);
        _map_centroid_to_polygon[centroid_cgal] = i; 
        _myocytes.push_back(std::move(myo));
        i++;
    }
    _tree = std::make_unique<Tree>(_points.begin(), _points.end());
}

substrate::substrate(const std::vector<std::string>& polygons, transform &t, std::function<double(double)> strain_f){

    _transform = std::move(t);
    _strain = strain_f;

    //iterate over polygons
    int i=0;
    auto start = std::chrono::high_resolution_clock::now(); // start timer

    for (const auto &poly: polygons){
        polyhedronSet myo(poly);
        Eigen::Vector3d centroid = myo.getUnstrainedCentroid();
        Kernel::Point_3 centroid_cgal(centroid(0), centroid(1), centroid(2));
        _points.emplace_back(centroid_cgal);
        _map_centroid_to_polygon[centroid_cgal] = i;
        _myocytes.push_back(std::move(myo));       
        i++;
    }
    
    _tree = std::make_unique<Tree>(_points.begin(), _points.end());
    auto end = std::chrono::high_resolution_clock::now(); // end timer
    
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time set substrate: " << elapsed.count() << " seconds" << std::endl;
}

void substrate::write_all_ply(const std::vector<polygon> &polyhedra, const std::string &filename)
{
    std::ofstream out(filename);
    if (!out.is_open())
    {
        throw std::runtime_error("Could not open file " + filename);
        return;
    }

    // Write header
    out << "ply\n";
    out << "format ascii 1.0\n";
    
    // Calculate total vertices and faces
    size_t totalVertices = 0;
    size_t totalFaces = 0;
    for (const auto &poly : polyhedra)
    {
        totalVertices += poly.getVertices().size();
        totalFaces += poly.getFaces().size();
    }

    out << "element vertex " << totalVertices << "\n";
    out << "property float x\n";
    out << "property float y\n";
    out << "property float z\n";
    out << "element face " << totalFaces << "\n";
    out << "property list uchar int vertex_indices\n";
    out << "end_header\n";

    // Write vertices
    for (const auto &poly : polyhedra)
    {
        for (auto &vertex : poly.getVertices())
        {
            out << vertex.x() << " " << vertex.y() << " " << vertex.z() << "\n";
        }
    }

    // Write faces
    size_t vertexOffset = 0;
    for (const auto &poly : polyhedra)
    {
        for (auto &face : poly.getFaces())
        {
            out << "3 " << face[0] + vertexOffset << " " << face[1] + vertexOffset << " " << face[2] + vertexOffset << "\n";
        }
        vertexOffset += poly.getVertices().size();
    }

    out.close();
}


void substrate::write_ply_surface_mesh(const Mesh &mesh, const std::string &filename)
{

    std::ofstream out(filename);
    if (!out.is_open())
    {
        throw std::runtime_error("Could not open file " + filename);
        return;
    }

    out << "ply\n";
    out << "format ascii 1.0\n";
    out << "element vertex " << mesh.number_of_vertices() << "\n";
    out << "property float x\n";
    out << "property float y\n";
    out << "property float z\n";
    out << "element face " << mesh.number_of_faces() << "\n";
    out << "property list uchar int vertex_indices\n";
    out << "end_header\n";

    // Write vertices
    for (const auto &v : mesh.vertices())
    {
        const auto &point = mesh.point(v);
        out << point.x() << " " << point.y() << " " << point.z() << "\n";
    }

    // Write faces
    for (const auto &f : mesh.faces())
    {
        out << "3";
        CGAL::Vertex_around_face_circulator<Mesh> vcirc(mesh.halfedge(f), mesh), done(vcirc);
        do
        {
            out << " " << (*vcirc).idx();
        } while (++vcirc != done);
        out << "\n";
    }

    out.close();
}

void substrate::write_ply_polyhedron(polygon &poly, const std::string &filename)
{
    std::ofstream out(filename);
    if (!out.is_open())
    {
        throw std::runtime_error("Could not open file " + filename);
        return;
    }

    // Write header
    out << "ply\n";
    out << "format ascii 1.0\n";
    out << "element vertex " << poly.getVertices().size() << "\n";
    out << "property float x\n";
    out << "property float y\n";
    out << "property float z\n";
    out << "element face " << poly.getFaces().size() << "\n";
    out << "property list uchar int vertex_indices\n";
    out << "end_header\n";

    // Write vertices
    for (auto &vertex : poly.getVertices())
    {
        out << vertex.x() << " " << vertex.y() << " " << vertex.z() << "\n";
    }

    for (auto &face : poly.getFaces())
    {
        out << "3 " << face[0] << " " << face[1] << " " << face[2] << "\n";
    }
    out.close();
}

boost::optional<std::tuple<int, double, Eigen::Vector3d>> substrate::intersectPolygon(const Eigen::Vector3d &point, const Eigen::Vector3d &step, int index_sequence) const
{
    // We will place the point in the center of the block as the tree is built with the centroid of the polygons

    // TODO: This will have to be modified when considering a full geometry. Maybe input into to the Ktree the 8 points of bounding box?
    // TODO: CHECK THAT THIS WORKS CURRENTLY, I DO NOT SEE WHY IT SHOULDNT!
    Kernel::Point_3 query((double)point(0), (double)point(1), 0.0);

    Neighbor_search search(*_trees[index_sequence], query, 30);

    double min_distance_children = std::numeric_limits<double>::max();
    int index_polygon_intersection = -1;
    Eigen::Vector3d normal_vector = Eigen::Vector3d::Zero();
    
    for (Neighbor_search::iterator it = search.begin(); it != search.end(); ++it)
    {
        //Pass in total_time here
    
        int index_polygon = _vector_map_centroid_to_polygon[index_sequence].at(it->first);
        const polygon& curr_poly = _myocytes[index_polygon].getPolygon(index_sequence);
        boost::optional<std::pair<int, double>> intersection = curr_poly.intersection(point, step);
        if (intersection)
        {
            int index_face = boost::get<std::pair<int, double>>(*intersection).first;
            Kernel::FT distance_to_face = boost::get<std::pair<int, double>>(*intersection).second;

            if (distance_to_face < min_distance_children)
            {
                min_distance_children = CGAL::to_double(distance_to_face);
                index_polygon_intersection = index_polygon;
                normal_vector = curr_poly.getNormalVector(index_face);
            }
        }
    }

    if (min_distance_children == std::numeric_limits<double>::max())
    {
        return boost::optional<std::tuple<int, double, Eigen::Vector3d>>();
    }
    
    return std::make_tuple(index_polygon_intersection, min_distance_children, normal_vector);
}

boost::optional<std::tuple<int, double, Eigen::Vector3d>>  substrate::intersectionBlock(const Eigen::Vector3d &point, const Eigen::Vector3d &step, int index_sequence) const
{
    
    // TODO: If we want reflecting boundaries in the block we need to output the normal vector of the face intersected
    CGAL::Bbox_3 solid_box = _transform.getBlockCurrentTime(index_sequence).getSolidBbox();
    CGAL::Point_3<Kernel> p1(point(0), point(1), point(2));
    CGAL::Point_3<Kernel> p2(point(0) + step(0), point(1) + step(1), point(2) + step(2));
    CGAL::Segment_3<Kernel> segment(p1, p2);

    if (p1.x() < solid_box.xmin() || p1.x() > solid_box.xmax() ||
        p1.y() < solid_box.ymin() || p1.y() > solid_box.ymax() ||
        p1.z() < solid_box.zmin() || p1.z() > solid_box.zmax())
    {
        //print eps between p1.x() and solid_box.xmin()
        //print eps between p1.y() and solid_box.ymin()
        //print eps between p1.z() and solid_box.zmin()
        //print eps between p1.x() and solid_box.xmax()
        //print eps between p1.y() and solid_box.ymax()
        //print eps between p1.z() and solid_box.zmax()
        std::cout << "eps between p1.x() and solid_box.xmin() = " << p1.x() - solid_box.xmin() << std::endl;
        std::cout << "eps between p1.y() and solid_box.ymin() = " << p1.y() - solid_box.ymin() << std::endl;
        std::cout << "eps between p1.z() and solid_box.zmin() = " << p1.z() - solid_box.zmin() << std::endl;
        std::cout << "eps between p1.x() and solid_box.xmax() = " << p1.x() - solid_box.xmax() << std::endl;
        std::cout << "eps between p1.y() and solid_box.ymax() = " << p1.y() - solid_box.ymax() << std::endl;
        std::cout << "eps between p1.z() and solid_box.zmax() = " << p1.z() - solid_box.zmax() << std::endl;

        throw std::runtime_error("how did we get here? p1 is outside the solid_box");
    }

    // Check if point p2 is inside the solid_box without the boundaries
    if (p2.x() > solid_box.xmin() && p2.x() < solid_box.xmax() &&
        p2.y() > solid_box.ymin() && p2.y() < solid_box.ymax() &&
        p2.z() > solid_box.zmin() && p2.z() < solid_box.zmax())
    {
        return boost::optional<std::tuple<int, double, Eigen::Vector3d>>();
    }

    boost::optional<std::pair<int, double>> intersection;
    try{
        intersection  = _transform.getBlockCurrentTime(index_sequence).intersection(point,step);
    }
    catch (std::exception& e){
        std::cout << "Substrate::intersectionBlock -> " << e.what() << std::endl;
        throw std::runtime_error("Substrate::intersectionBlock -> No intersection found but segment.target() is outside block");
    }

    if (intersection)
    {
            int index_face = boost::get<std::pair<int, double>>(*intersection).first;
            Kernel::FT distance_to_face = boost::get<std::pair<int, double>>(*intersection).second;
            Eigen::Vector3d normal_vector = _transform.getBlockCurrentTime(index_sequence).getNormalVector(index_face);
            double distance =  CGAL::to_double(distance_to_face);
            return std::make_tuple(index_face, distance, normal_vector);
    }
    else{
        throw std::runtime_error("WHATIS??");
    }   
}

int substrate::searchPolygon(const Eigen::Vector3d &point, int index_sequence, const std::string &frameOfReference ) const
{
    Eigen::Vector3d point_local;
    if (frameOfReference == "global")
    {
        double strain_z = strain_array_value(index_sequence);
        // TODO: Check if point is out of block and if it is check that point_local is. otherwise throw error means that transform is not init properly
        transform_info transform_data = getLocalFromGlobal(point,strain_z);
        point_local = transform_data.local_position;
    }
    else
    {
        point_local = point;
    }

    // TODO: Query the whole point when not using the block. Take into account when getting rid of block or implementing the transform type
    //  We will place the point in the center of the block as the tree is built with the centroid of the polygons
    Kernel::Point_3 query((double)point_local(0), (double)point_local(1), 0.0);

    Neighbor_search search(*_trees[index_sequence], query, 30);

    std::vector<std::pair<int, double>> indices_and_squared_distances;

    
    for (Neighbor_search::iterator it = search.begin(); it != search.end(); ++it)
    {
        int index_polygon = _vector_map_centroid_to_polygon[index_sequence].at(it->first);
        double squared_distance = CGAL::to_double(it->second); // Get the squared distance directly
        indices_and_squared_distances.push_back(std::make_pair(index_polygon, squared_distance));
    }

    
    std::sort(indices_and_squared_distances.begin(), indices_and_squared_distances.end(),
              [](const std::pair<int, double> &a, const std::pair<int, double> &b)
              {
                  return a.second < b.second;
              });

    // Check if the point is contained in any of the sorted polygons
    for (const auto &index_squared_distance_pair : indices_and_squared_distances)
    {
        int index_polygon = index_squared_distance_pair.first;
        const polygon& curr_poly = _myocytes[index_polygon].getPolygon(index_sequence);
        if (curr_poly.containsPoint(point_local))
        {
            return index_polygon;
        }
    }
    
    return -1;
}

void substrate::preComputeSubstrate(Eigen::VectorXd strain_array_time_dts)
{
    //Acumulate sum of sequence_dt;
    Eigen::VectorXd cumsum(strain_array_time_dts.size()+1);
    strain_array_value = Eigen::VectorXd(strain_array_time_dts.size()+1);
    cumsum(0) = 0;
    cumsum(1) = strain_array_time_dts(0);

    for (int i = 2; i < cumsum.size(); ++i) {
        cumsum(i) = cumsum(i - 1) + strain_array_time_dts(i-1);
    }

    strain_array_time = cumsum;

    //std::cout << "Strain array time: " << strain_array_time << std::endl;
    //Compute centroid of transform block

    for (int i = 0; i < strain_array_time.size(); ++i) {
        strain_array_value(i) = _strain(strain_array_time(i));
    }

    //std::cout << "Strain array value: " << strain_array_value << std::endl;

    Eigen::Vector3d centroid = _transform.get_block_centroid();


    //auto start = std::chrono::high_resolution_clock::now();
    //TODO: Move number of threads as input to the script , for now we hardcode it to 2
    #pragma omp parallel num_threads(1) shared(centroid,strain_array_time_dts,_strain)
    {
        #pragma omp master // This block will be executed by only one thread (the master)
        {
            std::cout << "Precomputing substrate for " << _myocytes.size() << " cells " << "with " << strain_array_time_dts.size() << " time steps" << std::endl;
            std::cout << "OpenMP is enabled. Number of threads: " << omp_get_num_threads() << std::endl;
            std::cout << "OpenMP version: " << _OPENMP << std::endl;
        }
        
        //TODO: static
        #pragma omp for schedule(static)
        for (int i = 0; i < _myocytes.size(); ++i) {
            _myocytes[i].precomputePolygons(centroid, strain_array_time_dts, _strain);
        }
    }

    //auto stop = std::chrono::high_resolution_clock::now();
    //auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    //std::cout << "Time taken by function: " << duration.count() << " milliseconds" << std::endl;



    /*for (auto &myo: _myocytes){
        myo.precomputePolygons(centroid,sequence_dt, _strain);
    }*/

    _vector_points.resize(strain_array_time.size());
    _vector_map_centroid_to_polygon.resize(strain_array_time.size());
    _trees.resize(strain_array_time.size());

    for (int i=0; i<strain_array_time.size(); ++i){

        for (int j=0; j<_myocytes.size(); ++j){
            Eigen::Vector3d centroid = _myocytes[j].getPolygon(i).getCentroid();
            Kernel::Point_3 centroid_cgal(centroid(0), centroid(1), 0.0);
            _vector_points[i].emplace_back(centroid_cgal);
            _vector_map_centroid_to_polygon[i][centroid_cgal] = j;
        }
        _trees[i] = std::make_unique<Tree>(_vector_points[i].begin(), _vector_points[i].end());
    }

    _transform.precomputeTransform(strain_array_time_dts, _strain);

    //Iterate over myocytes
    /*
    std::cout << "This is the length " <<  _myocytes[0]._polygonSet.size() << std::endl;
    for (int counter2 = 0; counter2 < _myocytes[0]._polygonSet.size(); ++counter2)
    {
        std::vector<polygon> all_polygons;  // Corrected declaration
        for (auto &myo : _myocytes)
        {
            std::cout << counter2 << std::endl;
            all_polygons.push_back(myo._polygonSet[counter2]);
        }
        std::string filename = "../../geometries_full/mesh_" + std::to_string(counter2) + ".ply";
        std::cout << filename << std::endl;
        //write_all_ply(all_polygons,filename);
    }*/
}

bool substrate::containsPoint(int index_polygon, const Eigen::Vector3d &point, int index_sequence)
{
    //return true;
    const polygon& curr_poly = _myocytes[index_polygon].getPolygon(index_sequence);
    return curr_poly.containsPoint(point);
}

Eigen::Vector3d substrate::getVoxelCentroid() const
{
    //From voxel get centroid
    Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
    centroid(0) = _voxel(0) + _voxel(1) / 2;
    centroid(1) = _voxel(2) + _voxel(3) / 2;
    centroid(2) = _voxel(4) + _voxel(5) / 2;
    return centroid;
}