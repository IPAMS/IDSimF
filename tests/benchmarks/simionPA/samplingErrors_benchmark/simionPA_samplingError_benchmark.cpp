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
#include "PSim_scalar_writer.hpp"
#include "appUtils_stopwatch.hpp"

#include <iostream>
//#include <string>

template<typename PAType>
void performBenchmark(unsigned int nSamples,
        double xStart,
        double xStop, double y, double z,
        std::string filename, std::string message){


    std::vector<double> xVec = ParticleSimulation::linspace(xStart,xStop, static_cast<int>(nSamples));

    PAType simPa(filename);

    ParticleSimulation::Scalar_writer scalarWriter("test.csv");

    std::cout << "Benchmark SIMION potential array " <<message << std::endl;
    AppUtils::Stopwatch stopWatch;
    stopWatch.start();

    for (unsigned int i=0; i< nSamples; ++i){

        //simPa.isElectrode(xVec[i],y,z);
        //simPa.getInterpolatedPotential(xVec[i],y,z);
        Core::Vector field = simPa.getField(xVec[i], y, z);
        double potential = simPa.getInterpolatedPotential(xVec[i], y, z);
        scalarWriter.writeTimestep(potential*10000.0, i);
        if (i % 10 == 0) {
            std::cout << "i:"<<i <<" x:"<< xVec[i]<< " field:" << field<<std::endl;
        }
    }

    stopWatch.stop();
    std::cout << "elapsed wall time:"<< stopWatch.elapsedSecondsWall()<<std::endl;
    std::cout << "elapsed cpu time:"<< stopWatch.elapsedSecondsCPU()<<std::endl;
}

int main() {

    Core::safetyGuards = true;

    //int n = 50000000;
    int n = 500;

    /*performBenchmark<ParticleSimulation::SimionPotentialArray>
            (n, 0.001, 0.038, 0.001, 0.002, "simion_test_planar_3d.pa","new");
    performBenchmark<ParticleSimulation::SimionPotentialArrayOriginal>
            (n, 0.001, 0.038, 0.001, 0.002, "simion_test_planar_3d.pa", "old");*/


    //performBenchmark<ParticleSimulation::SimionPotentialArrayOriginal>
    //    (n, 0.001, 10.38, 0.001, 0.002, "simion_test_cylindrical.pa", "old");
    performBenchmark<ParticleSimulation::SimionPotentialArray>
            (n, 0.0, 0.001, 0.001, 0.002, "simion_test_cylindrical.pa","new");

    return 0;
}
