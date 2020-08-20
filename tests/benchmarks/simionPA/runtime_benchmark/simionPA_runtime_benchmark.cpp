/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 2020 - Physical and Theoretical Chemistry /
 Institute of Pure and Applied Mass Spectrometry
 of the University of Wuppertal, Germany

 ------------

 A simple runtime benchmark for the SIMION Potential Array interface

 @author Walter Wissdorf
 ****************************/

#include "Core_debug.hpp"
#include "PSim_simionPotentialArray.hpp"
#include "PSim_math.hpp"

#include <ctime>
#include <iostream>
//#include <string>

template<typename PAType>
void performBenchmark(int nSamples,
        double xStart,
        double xStop, double y, double z,
        std::string filename, std::string message){


    std::vector<double> xVec = ParticleSimulation::linspace(xStart,xStop, nSamples);

    PAType simPa(filename);

    std::cout << "Benchmark SIMION potential array " <<message << std::endl;
    clock_t begin = std::clock();

    for (int i=0; i< nSamples; ++i){

        //simPa.isElectrode(xVec[i],y,z);
        //simPa.getInterpolatedPotential(xVec[i],y,z);
        simPa.getField(xVec[i], y, z);
    }

    clock_t end = std::clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "elapsed secs: "<< elapsed_secs<<std::endl;
}

int main(int argc, const char * argv[]) {

    Core::safetyGuards = true;

    int n = 50000000;

    performBenchmark<ParticleSimulation::SimionPotentialArray>
        (n, 0.038, 0.038, 0.001, 0.002, "simion_test_cylindrical_mirrored.pa","cylindrical");

    performBenchmark<ParticleSimulation::SimionPotentialArray>
            (n, -0.02, 0.02, 0.001, 0.0, "simion_test_planar_2d.pa","2d");

    performBenchmark<ParticleSimulation::SimionPotentialArray>
            (n, -0.036, 0.036, 0.001, 0.002, "simion_test_planar_3d_mirrored.pa","3d");

    return 0;
}
