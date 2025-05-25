cmake_minimum_required(VERSION 3.14)
project(NBodySimulation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(HDF5 REQUIRED COMPONENTS C CXX HL)

add_executable(nbody_cpp nbody_direct_sum_cpu.cpp)
target_include_directories(nbody_cpp PRIVATE ${HDF5_INCLUDE_DIRS})

if(UNIX AND NOT APPLE)
    # Sadece Linux'ta -pg ekle
    target_compile_options(nbody_cpp PRIVATE -pg)
    target_link_options(nbody_cpp PRIVATE -pg)
endif()

target_link_libraries(nbody_cpp PRIVATE ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES})