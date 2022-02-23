/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 2020 - Physical and Theoretical Chemistry /
 Institute of Pure and Applied Mass Spectrometry
 of the University of Wuppertal, Germany

 ------------

 A benchmark of the hard sphere collision model:
 A field free ion in a hard sphere gas bath should produce a maxwell boltzman
 velocity distribution

 @author Walter Wissdorf
 ****************************/

#include "CollisionModel_HardSphere.hpp"
#include "Core_particle.hpp"
#include "appUtils_stopwatch.hpp"

#include <ctime>
#include <iostream>

int main() {

    double diameterN2 = CollisionModel::HardSphereModel::DIAMETER_N2;
    CollisionModel::HardSphereModel hs = CollisionModel::HardSphereModel(10.0,298,28.0, diameterN2);
    Core::Particle ion = Core::Particle();

    ion.setDiameter(diameterN2);
    ion.setVelocity(Core::Vector(0,0,0));
    ion.setMassAMU(28.0);
    int n = 200000;

    std::cout << "Generate field free ion velocity samples" << std::endl;
    AppUtils::Stopwatch stopWatch;
    stopWatch.start();
    std::ofstream outputFile ("fieldfree_ion_samples.txt");

    for (int i =0; i<n; i++){
        hs.modifyVelocity(ion, 4e-7);
        outputFile << ion.getVelocity() << std::endl;
    }

    outputFile.close();

    stopWatch.stop();
    std::cout << "elapsed wall time:"<< stopWatch.elapsedSecondsWall()<<std::endl;
    std::cout << "elapsed cpu time:"<< stopWatch.elapsedSecondsCPU()<<std::endl;

    return 0;
}
