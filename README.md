## N-Body Simulation

An N-body simulation is a computational method used to study the dynamic evolution of a system of N interacting particles,
where each particle's motion is influenced by the forces from all other particles—typically under gravity. 
In this project, we implement a GPU-accelerated N-body simulation that models the positions, velocities, 
and gravitational interactions of multiple bodies in three-dimensional space.
The simulation calculates the gravitational forces between all particles at each time step, updates their positions 
and velocities accordingly, and handles inelastic collisions. Simulation results are periodically saved to HDF5 files 
for further analysis and visualization.

## Installation

To run gpu code:

```bash
sudo apt-get install libhdf5-dev
nvcc -o nbody_gpu nbody_direct_sum_gpu.cu -lhdf5_cpp -lhdf5
```