#include "PSim_trajectoryHDF5Writer.hpp"
#include "BTree_particle.hpp"
#include "appUtils_stopwatch.hpp"
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>



int main() {
    std::cout << "Benchmark HDF5 file writer" << std::endl;

    ParticleSimulation::TrajectoryHDF5Writer trajectoryWriter("test_hdf5_writer_trajectory.hd5");

    int nIons = 1e7;

    AppUtils::Stopwatch stopWatch;
    stopWatch.start();

    std::vector<std::unique_ptr<BTree::Particle>> particles;
    std::vector<BTree::Particle*> particlePtrs;

    for (int i=0; i<nIons; ++i){
        std::unique_ptr<BTree::Particle> newIon = std::make_unique<BTree::Particle>(
                Core::Vector(i*0.01, 0.0, 0.0), 1.0);

        newIon -> setSplatTime(i*0.1);
        particlePtrs.push_back(newIon.get());
        particles.push_back(std::move(newIon));
    }

    std::cout << "writing splattimes:: 33 "<<std::endl;
    trajectoryWriter.writeSplatTimes(particlePtrs);


    stopWatch.stop();

    std::cout << "elapsed wall time:"<< stopWatch.elapsedSecondsWall()<<std::endl;
    std::cout << "elapsed cpu time:"<< stopWatch.elapsedSecondsCPU()<<std::endl;
    return 0;
}
