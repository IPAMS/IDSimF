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

    clock_t begin = std::clock();

    for (int i =0; i<nSamples; ++i){
        hs.modifyVelocity(ion, 4e-6);
    }
    Core::Vector ionVelo = ion.getVelocity();


    clock_t end = std::clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    std::cout << "elapsed secs: "<< elapsed_secs<<std::endl;
    std::cout << "ion velocity: "<< ionVelo<<std::endl;
}

int main(int argc, const char * argv[]) {

    int n = 4000000;
    performBenchmark(n,false);
    performBenchmark(n,true);
    return 0;
}
