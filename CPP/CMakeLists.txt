cmake_minimum_required(VERSION 3.27)
project(NPSAT_URF)

#set(CMAKE_CXX_STANDARD 17)

#find_package(Eigen3 REQUIRED)

add_subdirectory(${DLIB_PATH} dlib_build)
include_directories(${EIGEN_PATH})

add_executable(NPSAT_URF main.cpp)

#target_link_libraries(NPSAT_URF Eigen3::Eigen) # dlib::dlib
target_link_libraries(NPSAT_URF dlib::dlib)