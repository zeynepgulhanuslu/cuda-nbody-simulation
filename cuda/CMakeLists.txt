cmake_minimum_required(VERSION 3.18)
project(nbody_gpu LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

find_package(HDF5 COMPONENTS CXX REQUIRED)
find_package(CUDAToolkit REQUIRED)

add_executable(nbody_gpu nbody_direct_sum_gpu.cu)

target_link_libraries(nbody_gpu
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
    CUDA::cudart
)

target_include_directories(nbody_gpu PRIVATE
    ${HDF5_INCLUDE_DIRS}
    ${CUDAToolkit_INCLUDE_DIRS}
)
