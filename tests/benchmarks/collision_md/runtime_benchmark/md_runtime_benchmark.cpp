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
#include "CLI11.hpp"

#include <omp.h>
#include <ctime>
#include <iostream>

void performBenchmark(size_t nSamples, size_t nParticles){

    double diameterHe = CollisionModel::MDInteractionsModel::DIAMETER_HE;
    FileIO::MolecularStructureReader reader = FileIO::MolecularStructureReader();
    std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection = reader.readMolecularStructure("test_molecularstructure_reader.csv");

    std::vector<Core::uniquePartPtr>particles;
    std::vector<Core::Particle*>particlesPtrs;
    double yPos = 0;
    for (size_t i=0; i<nParticles; ++i){
        Core::uniquePartPtr particle = std::make_unique<Core::Particle>(
                Core::Vector(0.0,yPos,0.0),
                Core::Vector(0.0,0.0,0.0),
                1.0,
                39);
        particle->setVelocity(Core::Vector(600.0, 0.0, 0.0));
        particle->setMolecularStructure(molecularStructureCollection.at("Ar+"));
        particle->setDiameter(particle->getMolecularStructure()->getDiameter());
        particlesPtrs.push_back(particle.get());
        particles.push_back(std::move(particle));
        yPos = yPos+0.0001;
    }
    // Core::Particle ion;
    // ion.setMolecularStructure(CollisionModel::MolecularStructure::molecularStructureCollection.at("H2+"));
    // ion.setVelocity(Core::Vector(100.0, 0.0, 0.0));
    // CollisionModel::MDInteractionsModel mdSim = CollisionModel::MDInteractionsModel(
    //                     10000,
    //                     298,
    //                     4.003,
    //                     diameterHe,
    //                     0.203e-30,
    //                     "He",
    //                     1e-10, 
    //                     1e-16);


    std::cout << "Benchmark molecular dynamics collision model "<<std::endl;

    AppUtils::Stopwatch stopWatch;
    stopWatch.start();

    std::size_t i;
    #pragma omp parallel \
        default(none) \
        firstprivate(particlesPtrs, nParticles, nSamples,  diameterHe, molecularStructureCollection)
    {
        CollisionModel::MDInteractionsModel mdSim = CollisionModel::MDInteractionsModel(
                10000,
                298,
                4.003,
                diameterHe,
                0.203e-30,
                "He",
                1e-10,
                1e-16,
                2,
                1,
                25, 
                molecularStructureCollection);

        for (size_t j = 0; j<nSamples; j++) {
            #pragma omp for schedule(dynamic, nParticles/20)
            for (i = 0; i<nParticles; ++i) {
                mdSim.modifyVelocity(*particlesPtrs[i], 1e-11);
            }
        }
    }
    stopWatch.stop();
    std::cout << "elapsed wall time:"<< stopWatch.elapsedSecondsWall()<<std::endl;
    std::cout << "elapsed cpu time:"<< stopWatch.elapsedSecondsCPU()<<std::endl;
}

int main(int argc, char** argv) {
    CLI::App app{"Simple benchmark of md collision model", "md collision model benchmark"};

    bool verbose = false;
    app.add_flag("-v,--verbose", verbose, "be verbose");

    unsigned int nParticles = 1;//23;
    app.add_option("--n_particles,-p", nParticles, "number of particles")->required();

    unsigned int nSamples = 1;
    app.add_option("--samples,-s", nSamples, "number of samples")->required();

    int numberOfThreads= 1;
    app.add_option("--n_threads,-t", numberOfThreads, "number of parallel threads")->required();
    CLI11_PARSE(app, argc, argv);

    omp_set_num_threads(numberOfThreads);
    performBenchmark(nSamples, nParticles);
    return 0;
}
