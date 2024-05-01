cmake_minimum_required(VERSION 3.27)
project(NPSAT_URF)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)

add_executable(NPSAT_URF main.cpp)

target_link_libraries(NPSAT_URF PRIVATE Eigen3::Eigen)