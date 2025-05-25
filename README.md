## N-Body Simulation

An N-body simulation is a computational method used to study the dynamic evolution of a system of N interacting particles,
where each particle's motion is influenced by the forces from all other particles—typically under gravity. 
In this project, we implement a GPU-accelerated N-body simulation that models the positions, velocities, 
and gravitational interactions of multiple bodies in three-dimensional space.
The simulation calculates the gravitational forces between all particles at each time step, updates their positions 
and velocities accordingly, and handles inelastic collisions. Simulation results are periodically saved to HDF5 files 
for further analysis and visualization.

First I implemented code using direct sum for CPU and GPU. I will improve my solution 
by implementing Barnes-Hut algorithm for GPU. 

## Installation
To run the CPU code:

```bash
cd cpu
mkdir build && cd build
cmake ..
make -j4
```
To run gpu code:

```bash
sudo apt-get install libhdf5-dev
mkdir build && cd build
cd build
cmake ..
make -j4
```

## Usage
To run the CPU simulation, execute the following command:

```bash
./nbody_cpu output_file number_of_bodies steps
```

To run the GPU simulation, execute the following command:

```bash
./nbody_gpu output_file number_of_bodies steps
```

## Profiling for CPU and GPU

To profile the CPU code, you can use `gprof`:'
```bash
gprof ./nbody_cpu gmon.out > cpu_profile.txt
```


To profile the GPU code, you can use `ncu`:
```bash
ncu --target-processes all ./nbody_gpu /media/N100_step300.h5 100 300
```


