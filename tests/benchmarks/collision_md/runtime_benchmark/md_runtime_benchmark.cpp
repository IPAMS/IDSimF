/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 2020 - Physical and Theoretical Chemistry /
 Institute of Pure and Applied Mass Spectrometry
 of the University of Wuppertal, Germany

 ------------

 A simple runtime benchmark for the Hard Sphere Collision Model in IDSimF

 @author Walter Wissdorf
 ****************************/

#include "CollisionModel_MDInteractions.hpp"
#include "FileIO_MolecularStructureReader.hpp"
#include "Core_particle.hpp"
#include "appUtils_stopwatch.hpp"

#include <ctime>
#include <iostream>

void performBenchmark(int nSamples){

    double diameterHe = CollisionModel::MDInteractionsModel::DIAMETER_HE;
    FileIO::MolecularStructureReader reader = FileIO::MolecularStructureReader();
    reader.readMolecularStructure("test_molecularstructure_reader.csv");
    Core::Particle ion;
    ion.setMolecularStructure(CollisionModel::MolecularStructure::molecularStructureCollection.at("H2+"));
    ion.setVelocity(Core::Vector(100.0, 0.0, 0.0));
    CollisionModel::MDInteractionsModel mdSim = CollisionModel::MDInteractionsModel(
                        100000,
                        298,
                        4.003,
                        2.89e-10,
                        0.203e-30,
                        "He",
                        500e-14, 
                        5e-14);


    std::cout << "Benchmark molecular dynamics collision model ";

    AppUtils::Stopwatch stopWatch;
    stopWatch.start();

    for (int i =0; i<nSamples; ++i){
        mdSim.modifyVelocity(ion, 4e-6);
    }
    Core::Vector ionVelo = ion.getVelocity();


    stopWatch.stop();
    std::cout << "elapsed wall time:"<< stopWatch.elapsedSecondsWall()<<std::endl;
    std::cout << "elapsed cpu time:"<< stopWatch.elapsedSecondsCPU()<<std::endl;
    std::cout << "ion velocity: "<< ionVelo<<std::endl;
}

int main() {

    int n = 4000;
    performBenchmark(n);
    return 0;
}
