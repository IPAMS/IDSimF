/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 2020 - Physical and Theoretical Chemistry /
 Institute of Pure and Applied Mass Spectrometry
 of the University of Wuppertal, Germany

 ------------

 A simple runtime benchmark for the Hard Sphere Collision Model in IDSimF

 @author Walter Wissdorf
 ****************************/

#include "CollisionModel_HardSphere.hpp"
#include "BTree_particle.hpp"
#include "appUtils_stopwatch.hpp"

#include <ctime>
#include <iostream>

void performBenchmark(int nSamples, bool maxwellApproximation){
    double diameterHe = CollisionModel::HardSphereModel::DIAMETER_HE;
    CollisionModel::HardSphereModel hs = CollisionModel::HardSphereModel(
            1.0,298,4.0,diameterHe,maxwellApproximation);
    BTree::Particle ion = BTree::Particle();
    ion.setDiameter(CollisionModel::HardSphereModel::DIAMETER_HE);
    ion.setVelocity(Core::Vector(100,0,0));
    ion.setMassAMU(28.0);


    std::cout << "Benchmark hard sphere collision model ";
    if (maxwellApproximation){
        std::cout << "with Maxwell approximation" << std::endl;
    }
    else{
        std::cout << "without Maxwell approximation" << std::endl;
    }

    AppUtils::Stopwatch stopWatch;
    stopWatch.start();

    for (int i =0; i<nSamples; ++i){
        hs.modifyVelocity(ion, 4e-6);
    }
    Core::Vector ionVelo = ion.getVelocity();


    stopWatch.stop();
    std::cout << "elapsed wall time:"<< stopWatch.elapsedSecondsWall()<<std::endl;
    std::cout << "elapsed cpu time:"<< stopWatch.elapsedSecondsCPU()<<std::endl;
    std::cout << "ion velocity: "<< ionVelo<<std::endl;
}

int main() {

    int n = 4000000;
    performBenchmark(n,false);
    performBenchmark(n,true);
    return 0;
}
