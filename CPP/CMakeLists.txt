cmake_minimum_required(VERSION 3.5)
project(NPSAT_URF LANGUAGES CXX)

find_package(MPI REQUIRED)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(EIGEN_PATH "" CACHE PATH "Path to the Eigen include directory")

if(NOT EIGEN_PATH)
    message(FATAL_ERROR "Eigen headers were not found. Configure with -DEIGEN_PATH=/path/to/eigen3")
endif()

file(TO_CMAKE_PATH "${EIGEN_PATH}" EIGEN_INCLUDE_DIR)

if(NOT EXISTS "${EIGEN_INCLUDE_DIR}/Eigen/Sparse")
    message(FATAL_ERROR "Eigen headers were not found under '${EIGEN_INCLUDE_DIR}'. Configure with -DEIGEN_PATH=/path/to/eigen3")
endif()

add_executable(NPSAT_URF main.cpp)
target_include_directories(NPSAT_URF PRIVATE
    "${EIGEN_INCLUDE_DIR}"
    "${MPI_CXX_INCLUDE_PATH}")
target_link_libraries(NPSAT_URF PRIVATE ${MPI_CXX_LIBRARIES})

if(MPI_CXX_COMPILE_FLAGS)
    set_target_properties(NPSAT_URF PROPERTIES
        COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
endif()

if(MPI_CXX_LINK_FLAGS)
    set_target_properties(NPSAT_URF PROPERTIES
        LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
endif()
