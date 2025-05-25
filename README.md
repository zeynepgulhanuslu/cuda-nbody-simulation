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
mkdir build
cd build
cmake ..
make -j4
```
To run gpu code:

```bash
sudo apt-get install libhdf5-dev
cd cuda
mkdir build
cd build
cmake ..
make -j4
```

## Usage
To run the CPU simulation, execute the following command:

```bash
./nbody_cpu -n <number_of_bodies> -t <simulation_time> -dt <time_step> -o <output_file>
```

To run the GPU simulation, execute the following command:

```bash
./nbody_gpu -n <number_of_bodies> -t <simulation_time> -dt <time_step> -o <output_file>
```


