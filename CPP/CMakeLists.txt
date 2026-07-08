cmake_minimum_required(VERSION 3.5)
project(NPSAT_URF)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(DEFINED EIGEN_PATH AND EXISTS "${EIGEN_PATH}/Eigen/Sparse")
    include_directories(${EIGEN_PATH})
elseif(EXISTS "/usr/include/eigen3/Eigen/Sparse")
    include_directories(/usr/include/eigen3)
else()
    message(FATAL_ERROR "Eigen headers were not found. Install Eigen or configure with -DEIGEN_PATH=/path/to/eigen3")
endif()

add_executable(NPSAT_URF main.cpp)
