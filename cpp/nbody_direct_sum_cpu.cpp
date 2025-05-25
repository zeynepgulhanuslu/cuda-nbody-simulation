#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <string>
#include "H5Cpp.h"
#include <sstream>
#include <filesystem>

namespace fs = std::filesystem;
using namespace H5;
using namespace std;

constexpr double G = 6.67430e-11;

struct Vec3 {
    double x, y, z;

    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }

    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }

    double norm() const { return std::sqrt(x*x + y*y + z*z); }
    Vec3 normalized() const { double n = norm(); return n > 0 ? (*this)/n : Vec3{0,0,0}; }
    double dot(const Vec3& o) const { return x*o.x + y*o.y + z*o.z; }
};

struct Body {
    Vec3 position;
    Vec3 velocity;
    Vec3 force;
    double mass;
};


void generateBodies(int N, vector<Body>& bodies,
                    double mass_min, double mass_max,
                    double pos_min, double pos_max,
                    double vel_min, double vel_max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> massDist(mass_min, mass_max);
    std::uniform_real_distribution<double> posDist(pos_min, pos_max);
    std::uniform_real_distribution<double> velDist(vel_min, vel_max);

    bodies.clear();
    bodies.reserve(N);

    for (int i = 0; i < N; i++) {
        Body b;
        b.mass = massDist(gen);
        b.position = {posDist(gen), posDist(gen), posDist(gen)};
        b.velocity = {velDist(gen), velDist(gen), velDist(gen)};
        b.force = {0,0,0};
        bodies.push_back(b);
    }
}

void computeForces(vector<Body>& bodies) {
    int N = bodies.size();

    for(auto& b : bodies) {
        b.force = {0,0,0};
    }

    for(int i=0; i<N; i++) {
        for(int j=0; j<N; j++) {
            if(i == j) continue;
            Vec3 diff = bodies[j].position - bodies[i].position;
            double dist = diff.norm() + 1e-10;
            double forceMag = G * bodies[i].mass * bodies[j].mass / (dist * dist);
            Vec3 forceVec = diff.normalized() * forceMag;
            bodies[i].force += forceVec;
        }
    }
}

void handleInelasticCollisions(vector<Body>& bodies, double radius, double e) {
    int N = bodies.size();
    for(int i=0; i<N; i++) {
        for(int j=i+1; j<N; j++) {
            Vec3 diff = bodies[j].position - bodies[i].position;
            double dist = diff.norm();
            if(dist < 2*radius) {
                Vec3 n_hat = diff.normalized();
                Vec3 v1 = bodies[i].velocity;
                Vec3 v2 = bodies[j].velocity;
                double m1 = bodies[i].mass;
                double m2 = bodies[j].mass;

                double v1_n = v1.dot(n_hat);
                double v2_n = v2.dot(n_hat);

                double v1_n_new = (m1 * v1_n + m2 * v2_n - m2 * e * (v1_n - v2_n)) / (m1 + m2);
                double v2_n_new = (m1 * v1_n + m2 * v2_n + m1 * e * (v1_n - v2_n)) / (m1 + m2);

                Vec3 v1_new = v1 + n_hat * (v1_n_new - v1_n);
                Vec3 v2_new = v2 + n_hat * (v2_n_new - v2_n);

                bodies[i].velocity = v1_new;
                bodies[j].velocity = v2_new;

                double overlap = 2*radius - dist;
                bodies[i].position -= n_hat * (overlap / 2);
                bodies[j].position += n_hat * (overlap / 2);
            }
        }
    }
}

void nbodyStep(vector<Body>& bodies, double dt, double radius, double e) {
    computeForces(bodies);
    for(auto& b : bodies) {
        Vec3 acceleration = b.force / b.mass;
        b.velocity += acceleration * dt;
        b.position += b.velocity * dt;
    }
    handleInelasticCollisions(bodies, radius, e);
}

bool writeAllStepsToHDF5(const std::string& filename, const std::vector<std::vector<Body>>& allSteps) {
    try {
        H5File file(filename, H5F_ACC_TRUNC);

        size_t numSteps = allSteps.size();
        size_t N = allSteps[0].size();
        Group stepGroup = file.createGroup("/steps");

        for (size_t step = 0; step < numSteps; step++) {
            std::ostringstream stepName;
            stepName << "step_" << std::setfill('0') << std::setw(3) << step;
            Group currentStep = stepGroup.createGroup(stepName.str());

            hsize_t dims[2] = {N, 3};

            // Positions dataset for this step
            DataSpace dataspace(2, dims);
            DataSet dataset_pos = currentStep.createDataSet("positions", PredType::NATIVE_DOUBLE, dataspace);
            std::vector<double> pos_data(N * 3);
            for (hsize_t i = 0; i < N; ++i) {
                pos_data[3 * i] = allSteps[step][i].position.x;
                pos_data[3 * i + 1] = allSteps[step][i].position.y;
                pos_data[3 * i + 2] = allSteps[step][i].position.z;
            }
            dataset_pos.write(pos_data.data(), PredType::NATIVE_DOUBLE);

            // Velocities dataset for this step
            DataSet dataset_vel = currentStep.createDataSet("velocities", PredType::NATIVE_DOUBLE, dataspace);
            std::vector<double> vel_data(N * 3);
            for (hsize_t i = 0; i < N; ++i) {
                vel_data[3 * i] = allSteps[step][i].velocity.x;
                vel_data[3 * i + 1] = allSteps[step][i].velocity.y;
                vel_data[3 * i + 2] = allSteps[step][i].velocity.z;
            }
            dataset_vel.write(vel_data.data(), PredType::NATIVE_DOUBLE);

            // Masses dataset for this step (1D)
            hsize_t dims_mass[1] = {N};
            DataSpace mass_space(1, dims_mass);
            DataSet dataset_mass = currentStep.createDataSet("masses", PredType::NATIVE_DOUBLE, mass_space);
            std::vector<double> mass_data(N);
            for (hsize_t i = 0; i < N; ++i) {
                mass_data[i] = allSteps[step][i].mass;
            }
            dataset_mass.write(mass_data.data(), PredType::NATIVE_DOUBLE);
        }

        // Save simulation info
        Group infoGroup = file.createGroup("/info");

        hsize_t scalar_dims[1] = {1};
        DataSpace scalar_space(1, scalar_dims);

        DataSet dataset_numSteps = infoGroup.createDataSet("num_steps", PredType::NATIVE_INT, scalar_space);
        int numStepsInt = static_cast<int>(numSteps);
        dataset_numSteps.write(&numStepsInt, PredType::NATIVE_INT);

        DataSet dataset_numParticles = infoGroup.createDataSet("num_particles", PredType::NATIVE_INT, scalar_space);
        int numParticlesInt = static_cast<int>(N);
        dataset_numParticles.write(&numParticlesInt, PredType::NATIVE_INT);

        return true;
    } catch (FileIException& e) {
        e.printErrorStack();
        return false;
    } catch (DataSetIException& e) {
        e.printErrorStack();
        return false;
    } catch (DataSpaceIException& e) {
        e.printErrorStack();
        return false;
    }
}

int main(int argc, char* argv[]) {
    int N = 15;
    int steps = 300;
    double dt = 1e6;
    double radius = 2e9;
    double e = 0.7;

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " output_file [N steps]" << std::endl;
        return 1;
    }
    std::string output_file = argv[1];
    if (argc > 2) N = std::atoi(argv[2]);
    if (argc > 3) steps = std::atoi(argv[3]);

    // create directory if needed
    fs::path out_path(output_file);
    if (!out_path.parent_path().empty()) {
        try {
            fs::create_directories(out_path.parent_path());
        } catch (const fs::filesystem_error& e) {
            std::cerr << "Error occur when creating output directory: " << e.what() << std::endl;
            return 1;
        }
    }

    std::vector<Body> bodies;
    generateBodies(N, bodies, 1e20, 1e25, -1e11, 1e11, -1e3, 1e3);
    std::vector<std::vector<Body>> allSteps;
    allSteps.push_back(bodies);
    auto start = std::chrono::high_resolution_clock::now();

    for (int step = 0; step < steps; ++step) {
        nbodyStep(bodies, dt, radius, e);
        allSteps.push_back(bodies);
        if ((step+1) % 10 == 0 || step == steps-1) {
            std::cout << "Progress: " << step+1 << "/" << steps << " (%"
                      << std::fixed << std::setprecision(1) << (step+1)*100.0/steps << ")\n";
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "All steps are calculated. Saving to HDF5 file.\n";

    if (!writeAllStepsToHDF5(output_file, allSteps)) {
        std::cerr << "An error occurred when saving results.\n";
        return 1;
    }

    std::cout << "========\n";
    std::cout << "Process completed! Total time: " << elapsed_seconds.count() << " seconds\n";
    std::cout << "Results saved to " << output_file << " file\n";
    std::cout << "========\n";

    return 0;
}
