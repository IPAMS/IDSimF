/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 2020 - Physical and Theoretical Chemistry /
 Institute of Pure and Applied Mass Spectrometry
 of the University of Wuppertal, Germany

 ------------

 A simple runtime benchmark for the Statistical Diffusion Model in IDSimF

 @author Walter Wissdorf
 ****************************/

#include "CollisionModel_StatisticalDiffusion.hpp"
#include "BTree_particle.hpp"

#include <ctime>
#include <iostream>

void performBenchmark(int nSamples, double dt){

    std::cout << "Benchmark SDS collision model dt="<<dt << std::endl;
    clock_t begin = std::clock();

    CollisionModel::StatisticalDiffusionModel sds(100000,298,28,0.366 * 1.0e-9);
    BTree::Particle ion;
    ion.setMassAMU(100);
    ion.setVelocity({100,0.0,0.0});
    sds.setSTPParameters(ion);
    sds.updateModelParameters(ion);

    for (int i=0; i< nSamples; ++i){
        Core::Vector acceleration {200, 0.0, 0.0};
        sds.modifyAcceleration(acceleration, ion, dt);
        sds.modifyPosition(ion.getLocation(), ion, dt);
    }

    Core::Vector ionLoc = ion.getLocation();

    clock_t end = std::clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    std::cout << "ion location: "<< ionLoc << std::endl;
    std::cout << "elapsed secs: "<< elapsed_secs<<std::endl;
}

int main(int argc, const char * argv[]) {

    int n = 4000000;
    performBenchmark(n, 1e-2);
    performBenchmark(n, 1e-1);
    return 0;
}
