cmake_minimum_required(VERSION 3.5)
set(CMAKE_C_COMPILER "/usr/bin/gcc")
set(CMAKE_CXX_COMPILER "/usr/bin/g++")
set(CMAKE_CXX_STANDARD 14)
project(3DRandomWalk)

# If you need to debug, set CMAKE_BUILD_TYPE to "Debug" and optimize to "-O0" and run locally without Docker. You will have to install dependencies 
# through anaconda or directly in the system.
set(CMAKE_BUILD_TYPE "Release")
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # Choose your desired optimization level and floating-point settings
    set(OPTIMIZATION_LEVEL "-O3")
    set(FP_CONSISTENCY_FLAGS "")

    # Apply the same optimization level and floating-point settings to Debug and Release configurations
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OPTIMIZATION_LEVEL} ${FP_CONSISTENCY_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${OPTIMIZATION_LEVEL} ${FP_CONSISTENCY_FLAGS}")
endif()


#IF YOU USE ANACONDA ENV uncomment the following lines:
    # set(CONDA_ENV_NAME "randomwalk")

    # # Set the paths to the conda environment
    # set(CONDA_PREFIX "~/anaconda3/envs/${CONDA_ENV_NAME}")
    # set(CONDA_INCLUDE_DIR "${CONDA_PREFIX}/include")
    # set(CONDA_LIBRARY_DIR "${CONDA_PREFIX}/lib")

    # #set(Eigen3_DIR "${CONDA_PREFIX}/share/eigen3/cmake")
    # set(CMAKE_PREFIX_PATH ${CONDA_PREFIX} ${CMAKE_PREFIX_PATH})


set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${CMAKE_CURRENT_SOURCE_DIR}/third_party/yaml-cpp")

find_package(Eigen3 3.3 REQUIRED)
message(STATUS "Eigen3 found at: ${Eigen3_DIR}")

find_package(CGAL REQUIRED COMPONENTS Core)
message(STATUS "cgal found at: ${CGAL_DIR}")

find_package(Boost REQUIRED COMPONENTS filesystem system)
message(STATUS "Boost found at: ${Boost_DIR}")

find_package(Catch2 REQUIRED)
message(STATUS "Catch2 found at: ${Catch2_DIR}")

find_package(matioCpp REQUIRED)
message(STATUS "matioCpp found at: ${matioCpp_DIR}")

find_package(OpenMP)
message(STATUS "OpenMP found at: ${OpenMP_CXX_FLAGS}")

find_package(Ceres REQUIRED)
message(STATUS "Ceres found at: ${Ceres_DIR}")

set(GEOMETRY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Geometry")
set(SUBSTRATE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Substrate")
set(MONTECARLO_DIR "${CMAKE_CURRENT_SOURCE_DIR}/MonteCarlo")
set(MRI_DIR "${CMAKE_CURRENT_SOURCE_DIR}/MRI")
set(UTILITY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Utility")

include_directories(${CERES_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/third_party/yaml-cpp/include)

file(GLOB_RECURSE SOURCES "${GEOMETRY_DIR}/*.cpp" "${SUBSTRATE_DIR}/*.cpp" "${MONTECARLO_DIR}/*.cpp" "${MRI_DIR}/*.cpp" "${UTILITY_DIR}/*.cpp" "SimulationRunner.cpp")
file(GLOB_RECURSE HEADERS "${GEOMETRY_DIR}/*.h" "${SUBSTRATE_DIR}/*.h" "${MONTECARLO_DIR}/*.h" "${MRI_DIR}/*.h" "${UTILITY_DIR}/*.h" "SimulationRunner.h")

# Add the source files for your project
#add_executable(myTests testing/oneGenerator.h testing/testing_simulation.cpp testing/testing_sequence.cpp testing/testing_polygon.cpp testing/testing_substrate.cpp ${SOURCES} ${HEADERS})
add_executable(3DRandomWalk run_sim.cpp ${SOURCES} ${HEADERS})

# Setting RPATH for 3DRandomWalk target
# Link against the yaml-cpp library
target_link_libraries(3DRandomWalk Eigen3::Eigen Catch2::Catch2WithMain OpenMP::OpenMP_CXX ${Boost_LIBRARIES} CGAL::CGAL matioCpp::matioCpp Boost::filesystem Boost::system yaml-cpp ${CERES_LIBRARIES})
target_compile_definitions(3DRandomWalk PRIVATE CGAL_PMP_USE_CERES_SOLVER )