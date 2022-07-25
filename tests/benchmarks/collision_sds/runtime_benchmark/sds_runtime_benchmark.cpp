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
#include "Core_particle.hpp"
#include "appUtils_stopwatch.hpp"

#include <iostream>

void performBenchmark(int nSamples, double dt){

    std::cout << "Benchmark SDS collision model dt="<<dt << std::endl;
    AppUtils::Stopwatch stopWatch;
    stopWatch.start();

    CollisionModel::StatisticalDiffusionModel sds(100000,298,28,0.366 * 1.0e-9);
    Core::Particle ion;
    ion.setMassAMU(100);
    ion.setVelocity({100,0.0,0.0});
    sds.setSTPParameters(ion);
    sds.updateModelParticleParameters(ion);

    for (int i=0; i< nSamples; ++i){
        Core::Vector acceleration {200, 0.0, 0.0};
        sds.modifyAcceleration(acceleration, ion, dt);
        sds.modifyPosition(ion.getLocation(), ion, dt);
    }

    Core::Vector ionLoc = ion.getLocation();

    stopWatch.stop();

    std::cout << "ion location: "<< ionLoc << std::endl;
    std::cout << "elapsed wall time:"<< stopWatch.elapsedSecondsWall()<<std::endl;
    std::cout << "elapsed cpu time:"<< stopWatch.elapsedSecondsCPU()<<std::endl;

}

int main() {

    int n = 4000000;
    performBenchmark(n, 1e-2);
    performBenchmark(n, 1e-1);
    return 0;
}
