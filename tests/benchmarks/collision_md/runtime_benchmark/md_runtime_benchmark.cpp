/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 2020 - Physical and Theoretical Chemistry /
 Institute of Pure and Applied Mass Spectrometry
 of the University of Wuppertal, Germany

 ------------

 A simple runtime benchmark for the MD Collision Model
 ****************************/

#include "CollisionModel_MDInteractions.hpp"
#include "CollisionModel_HardSphere.hpp"
#include "FileIO_MolecularStructureReader.hpp"
#include "Core_particle.hpp"
#include "appUtils_stopwatch.hpp"
#include "CLI11.hpp"

#include <omp.h>
#include <ctime>
#include <iostream>

// Implement simple test method
//void

void myWorkFct(double dt) {
    //Core::RandomSource* rndSource = Core::globalRandomGeneratorPool->getThreadRandomSource();

    // Calculate collision cross section between particle and collision gas:
    //   TODO: It seems to be unnecessary to constantly recalculate this
    //   value, cache the calculated values somehow?
    Core::Vector velocity = {0.1,0.1,0.1};

    double x = velocity.x();
    double y = velocity.y();

    for (int i=0; i<100; ++i){
        x = x * sin(i*dt);
        y = y * cos(i*dt);
    }

    velocity.x(x);
    velocity.y(y);

//ion.setVelocity(velocity);
}


void performBenchmark(size_t nSamples, size_t nParticles){

    //double diameterHe = CollisionModel::MDInteractionsModel::DIAMETER_HE;
    //FileIO::MolecularStructureReader reader = FileIO::MolecularStructureReader();
    //reader.readMolecularStructure("test_molecularstructure_reader.csv");

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
        //particle->setMolecularStructure(CollisionModel::MolecularStructure::molecularStructureCollection.at("Ar+"));
        //particle->setDiameter(particle->getMolecularStructure()->getDiameter());
        particlesPtrs.push_back(particle.get());
        particles.push_back(std::move(particle));
        yPos = yPos+0.0001;
    }
    // Core::Particle ion;
    // ion.setMolecularStructure(CollisionModel::MolecularStructure::molecularStructureCollection.at("H2+"));
    // ion.setVelocity(Core::Vector(100.0, 0.0, 0.0));



    std::cout << "Benchmark molecular dynamics collision model  << Iter 4"<<std::endl;

    AppUtils::Stopwatch stopWatch;
    stopWatch.start();


    //shared(diameterHe, particlesPtrs)
    //{
        //Core::Particle myParticle = *particlesPtrs[0];
        /*CollisionModel::MDInteractionsModel mdSim = CollisionModel::MDInteractionsModel(
                10000,
                298,
                4.003,
                diameterHe,
                0.203e-30,
                "He",
                1e-10,
                1e-16);*/

        /*CollisionModel::HardSphereModel hs = CollisionModel::HardSphereModel(
                1.0,298,4.0,diameterHe,true);*/

        for (size_t j = 0; j<nSamples; j++) {
            std::size_t i;
            #pragma omp parallel \
                default(none) \
                private(i) \
                firstprivate(nParticles, nSamples)
                {
                    #pragma omp for
                    for (i = 0; i<nParticles; ++i) {
                        //mdSim.modifyVelocity(*particlesPtrs[i], 1e-11);
                        //hs.modifyVelocity(*particlesPtrs[i], 1e-11);
                        myWorkFct(1e-11);
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
    performBenchmark(nParticles, nSamples);
    return 0;
}
