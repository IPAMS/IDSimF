#include "BTree_tree.hpp"
#include "Core_particle.hpp"
#include "Integration_verletIntegrator.hpp"
#include "Integration_parallelVerletIntegrator.hpp"
#include "Integration_fmm3dIntegrator.hpp"
#include "FileIO_trajectoryHDF5Writer.hpp"
#include "PSim_util.hpp"
#include "CollisionModel_StatisticalDiffusion.hpp"
#include "appUtils_stopwatch.hpp"
#include "CLI11.hpp"
#include <iostream>
#include <numeric>
#include <omp.h>


void runIntegrator(Integration::AbstractTimeIntegrator &integrator, unsigned int timeSteps, unsigned int numberOfParticles, double dt, std::string message){
    std::cout << "Benchmark " <<message << std::endl;
    std::cout << "number of particles: " <<numberOfParticles << std::endl;
    AppUtils::Stopwatch stopWatch;
    stopWatch.start();

    integrator.run(timeSteps, dt);

    stopWatch.stop();
    std::cout << "elapsed wall time:"<< stopWatch.elapsedSecondsWall()<<std::endl;
    std::cout << "elapsed cpu time:"<< stopWatch.elapsedSecondsCPU()<<std::endl;

}


unsigned int prepareIons(std::vector<std::unique_ptr<Core::Particle>> &particles,
                 std::vector<Core::Particle*> &particlePtrs, unsigned int nPerDirection){

    unsigned int nTotal = 0;
    for (unsigned int i=0; i<nPerDirection; i++){
        double pX = i*1.0/nPerDirection;
        for (unsigned int j=0; j<nPerDirection; j++){
            double pY = j*1.0/nPerDirection;
            for (unsigned int k=0; k<nPerDirection; k++){
                double pZ = k*1.0/nPerDirection;
                std::unique_ptr<Core::Particle> newIon = std::make_unique<Core::Particle>(Core::Vector(pX,pY,pZ), 1.0);
                newIon -> setMassAMU(100);
                particlePtrs.push_back(newIon.get());
                particles.push_back(std::move(newIon));
                nTotal++;
            }
        }
    }
    return nTotal;
}

template<class integratorT> std::vector<Core::Vector> runSimulation(unsigned int nIonsPerDirection, unsigned int timeSteps,
                                                                    double dt, double spaceChargeFactor,
                                                                    bool useCollisionModel, std::string runName){

    // define functions for the trajectory integration ==================================================
    auto accelerationFunction =
            [spaceChargeFactor](
                    Core::Particle *particle, int /*particleIndex*/,
                    SpaceCharge::FieldCalculator& scFieldCalculator, double /*time*/, int /*timestep*/) -> Core::Vector{

                double particleCharge = particle->getCharge();

                Core::Vector spaceChargeForce(0,0,0);
                if (spaceChargeFactor > 0) {
                    Core::Vector sf = scFieldCalculator.getEFieldFromSpaceCharge(*particle);
                    spaceChargeForce =
                            sf * (particleCharge * spaceChargeFactor);
                }

                return (spaceChargeForce / particle->getMass());
            };

    std::vector<std::unique_ptr<Core::Particle>> particles;
    std::vector<Core::Particle*>particlePtrs;

    unsigned int nIonsTotal = prepareIons(particles, particlePtrs, nIonsPerDirection);

    CollisionModel::StatisticalDiffusionModel sdsCollisionModel(100000.0, 298, 28, 3.64e-9);
    CollisionModel::AbstractCollisionModel *collisionModel;
    if (useCollisionModel){
        collisionModel = &sdsCollisionModel;
        for (unsigned int i=0; i<nIonsTotal; ++i){
            sdsCollisionModel.setSTPParameters(*particlePtrs[i]);
        }
    }
    else {
        collisionModel = nullptr;
    }

    integratorT integrator(
            particlePtrs,
            accelerationFunction, nullptr, nullptr, nullptr, collisionModel);

    runIntegrator(integrator, timeSteps, nIonsTotal, dt, runName);

    std::vector<Core::Vector> result;
    for (unsigned int i=0; i<nIonsTotal; ++i){
        result.push_back(particles[i]->getLocation());
    }

    return result;
}

int main(int argc, char** argv) {
    CLI::App app{"Simple benchmark of space charge calculation", "simpleSpaceCharge benchmark"};

    //app.add_option("-f,--file", filename, "A help string");

    bool useCollisionModel = false;
    app.add_flag("-c,--collisionModel", useCollisionModel, "Use collision model");

    bool verbose = false;
    app.add_flag("-v,--verbose", verbose, "be verbose");


    unsigned int nIonsPerDirection = 1;//23;
    app.add_option("--n_ions,-i", nIonsPerDirection, "number of ions per direction")->required();

    unsigned int timeSteps = 1;
    app.add_option("--time_steps,-t", timeSteps, "number of time steps")->required();
    double dt = 1e-3;
    double spaceChargeFactor = 1.0;

    int numberOfThreads_= 1;
    app.add_option("--n_threads,-n", numberOfThreads_, "number of parallel threads")->required();
    CLI11_PARSE(app, argc, argv);

    omp_set_num_threads(numberOfThreads_);

    //auto hdf5Writer = std::make_unique<FileIO::TrajectoryHDF5Writer>(
    //        "test_trajectories.hd5");


    // simulate ===============================================================================================
    //std::vector<Core::Vector> locationsSerial = runSimulation<Integration::VerletIntegrator>(
    //        nIonsPerDirection, timeSteps, dt, spaceChargeFactor, useCollisionModel, "serial");

    std::vector<Core::Vector> locationsParallel = runSimulation<Integration::ParallelVerletIntegrator>(
            nIonsPerDirection, timeSteps, dt, spaceChargeFactor, useCollisionModel, "parallel BTree");

    std::vector<Core::Vector> locationsFmm = runSimulation<Integration::FMMVerletIntegrator>(
            nIonsPerDirection, timeSteps, dt, spaceChargeFactor, useCollisionModel, "FMM3D");

    /*std::size_t nIonsTotal = locationsParallel.size();

    //std::vector<double> diffMagsTrees;
    std::vector<double> diffMagsFMM;
    for (unsigned int i=0; i<nIonsTotal; ++i){
        //diffMagsTrees.push_back( (locationsSerial[i] - locationsParallel[i]).magnitude() );
        diffMagsFMM.push_back( (locationsParallel[i] - locationsFmm[i]).magnitude() );
    }
    //double sumTrees = std::accumulate(diffMagsTrees.begin(), diffMagsTrees.end(), 0.0);
    double sumFmm = std::accumulate(diffMagsFMM.begin(), diffMagsFMM.end(), 0.0);

    if (verbose) {
        for (unsigned int i = 0; i<nIonsTotal; ++i) {
            std::cout << locationsParallel[i] << " | " << locationsFmm[i] << " | "
                      << (locationsParallel[i]-locationsFmm[i]).magnitude()
                      << std::endl;
        }
    }
    //std::cout << "sum diff trees: " << sumTrees << "sum diff fmm: " << sumFmm <<std::endl;
    std::cout << "sum diff fmm: " << sumFmm <<std::endl;*/
}