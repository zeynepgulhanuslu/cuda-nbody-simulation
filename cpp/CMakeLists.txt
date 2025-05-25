cmake_minimum_required(VERSION 3.14)
project(NBodySimulation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(HDF5 REQUIRED COMPONENTS C CXX HL)
set(USE_EXPERIMENTAL_FILESYSTEM OFF)
add_executable(nbody nbody_direct_sum_cpu.cpp) 

target_include_directories(nbody PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(nbody PRIVATE ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES})

