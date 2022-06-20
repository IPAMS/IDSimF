/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 2020 - Physical and Theoretical Chemistry /
 Institute of Pure and Applied Mass Spectrometry
 of the University of Wuppertal, Germany

 ------------

 A simple runtime benchmark for the MD Collision Model
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

    std::vector<Core::uniquePartPtr>particles;
    std::vector<Core::Particle*>particlesPtrs;
    size_t nParticles = 200;
    double yPos = 0;
    for (size_t i=0; i<nParticles; ++i){
        Core::uniquePartPtr particle = std::make_unique<Core::Particle>(
                Core::Vector(0.0,yPos,0.0),
                Core::Vector(0.0,0.0,0.0),
                1.0,
                39);
        particle->setVelocity(Core::Vector(600.0, 0.0, 0.0));
        particle->setMolecularStructure(CollisionModel::MolecularStructure::molecularStructureCollection.at("Ar+"));
        particle->setDiameter(particle->getMolecularStructure()->getDiameter());
        particlesPtrs.push_back(particle.get());
        particles.push_back(std::move(particle));
        yPos = yPos+0.0001;
    }
    // Core::Particle ion;
    // ion.setMolecularStructure(CollisionModel::MolecularStructure::molecularStructureCollection.at("H2+"));
    // ion.setVelocity(Core::Vector(100.0, 0.0, 0.0));
    CollisionModel::MDInteractionsModel mdSim = CollisionModel::MDInteractionsModel(
                        10000,
                        298,
                        4.003,
                        diameterHe,
                        0.203e-30,
                        "He",
                        1e-10, 
                        1e-16);


    std::cout << "Benchmark molecular dynamics collision model ";

    AppUtils::Stopwatch stopWatch;
    stopWatch.start();

    for(int j = 0; j < nSamples; j++){
        #pragma omp parallel for num_threads(7)
        for (size_t i =0; i<nParticles; ++i){
            mdSim.modifyVelocity(*particlesPtrs[i], 1e-11);
        }
    }


    stopWatch.stop();
    std::cout << "elapsed wall time:"<< stopWatch.elapsedSecondsWall()<<std::endl;
    std::cout << "elapsed cpu time:"<< stopWatch.elapsedSecondsCPU()<<std::endl;
  
}

int main() {

    int n = 400;
    performBenchmark(n);
    return 0;
}
