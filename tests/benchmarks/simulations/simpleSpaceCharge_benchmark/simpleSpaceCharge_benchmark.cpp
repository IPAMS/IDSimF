#include "BTree_tree.hpp"
#include "BTree_parallelTree.hpp"
#include "BTree_particle.hpp"
#include "Integration_verletIntegrator.hpp"
#include "Integration_parallelVerletIntegrator.hpp"
#include "PSim_trajectoryHDF5Writer.hpp"
#include "PSim_util.hpp"
#include "CollisionModel_StatisticalDiffusion.hpp"
#include "appUtils_stopwatch.hpp"
#include "CLI11.hpp"
#include <iostream>
#include <numeric>

void runIntegrator(Integration::AbstractTimeIntegrator &integrator, unsigned int timeSteps, double dt, std::string message){
    std::cout << "Benchmark " <<message << std::endl;
    AppUtils::Stopwatch stopWatch;
    stopWatch.start();

    integrator.run(timeSteps, dt);

    stopWatch.stop();
    std::cout << "elapsed wall time:"<< stopWatch.elapsedSecondsWall()<<std::endl;
    std::cout << "elapsed cpu time:"<< stopWatch.elapsedSecondsCPU()<<std::endl;

}


unsigned int prepareIons(std::vector<std::unique_ptr<BTree::Particle>> &particles,
                 std::vector<BTree::Particle*> &particlePtrs, unsigned int nPerDirection){

    unsigned int nTotal = 0;
    for (unsigned int i=0; i<nPerDirection; i++){
        double pX = i*1.0/nPerDirection;
        for (unsigned int j=0; j<nPerDirection; j++){
            double pY = j*1.0/nPerDirection;
            for (unsigned int k=0; k<nPerDirection; k++){
                double pZ = k*1.0/nPerDirection;
                std::unique_ptr<BTree::Particle> newIon = std::make_unique<BTree::Particle>(Core::Vector(pX,pY,pZ), 1.0);
                newIon -> setMassAMU(100);
                particlePtrs.push_back(newIon.get());
                particles.push_back(std::move(newIon));
                nTotal++;
            }
        }
    }
    return nTotal;
}

int main(int argc, char** argv) {
    CLI::App app{"Simple benchmark of space charge calculation", "simpleSpaceCharge benchmark"};

    //app.add_option("-f,--file", filename, "A help string");

    bool useCollisionModel = false;
    app.add_flag("-c,--collisionModel", useCollisionModel, "Use collision model");

    bool verbose = false;
    app.add_flag("-v,--verbose", verbose, "be verbose");
    CLI11_PARSE(app, argc, argv);


    unsigned int nIonsPerDirection = 23;
    unsigned int timeSteps = 200;
    double dt = 1e-3;
    double spaceChargeFactor = 1.0;

    // define functions for the trajectory integration ==================================================
    auto accelerationFunctionSerial =
            [spaceChargeFactor](
                    BTree::Particle *particle, int /*particleIndex*/,
                    BTree::Tree &tree, double /*time*/, int /*timestep*/) -> Core::Vector{

                double particleCharge = particle->getCharge();

                Core::Vector spaceChargeForce(0,0,0);
                if (spaceChargeFactor > 0) {
                    spaceChargeForce =
                            tree.computeEFieldFromTree(*particle) * (particleCharge * spaceChargeFactor);
                }
                return (spaceChargeForce / particle->getMass());
            };


    auto accelerationFunctionParallelNew =
            [spaceChargeFactor](
                    BTree::Particle *particle, int /*particleIndex*/,
                    BTree::ParallelTree &tree, double /*time*/, int /*timestep*/) -> Core::Vector{

                double particleCharge = particle->getCharge();

                Core::Vector spaceChargeForce(0,0,0);
                if (spaceChargeFactor > 0) {
                    spaceChargeForce =
                            tree.computeEFieldFromTree(*particle) * (particleCharge * spaceChargeFactor);
                }
                return (spaceChargeForce / particle->getMass());
            };

    auto hdf5Writer = std::make_unique<FileIO::TrajectoryHDF5Writer>(
            "test_trajectories.hd5");

    std::vector<std::unique_ptr<BTree::Particle>> particlesSerial;
    std::vector<BTree::Particle*>particlePtrsSerial;
    std::vector<std::unique_ptr<BTree::Particle>> particlesParallelNew;
    std::vector<BTree::Particle*>particlePtrsParallelNew;

    unsigned int nIonsTotal = prepareIons(particlesSerial, particlePtrsSerial, nIonsPerDirection);
    prepareIons(particlesParallelNew, particlePtrsParallelNew, nIonsPerDirection);
    
    CollisionModel::StatisticalDiffusionModel sdsCollisionModel(100000.0, 298, 28, 3.64e-9);
    CollisionModel::AbstractCollisionModel *collisionModel;
    if (useCollisionModel){
        collisionModel = &sdsCollisionModel;
        for (unsigned int i=0; i<nIonsTotal; ++i){
            sdsCollisionModel.setSTPParameters(*particlePtrsSerial[i]);
            sdsCollisionModel.setSTPParameters(*particlePtrsParallelNew[i]);
        }

    }
    else {
        collisionModel = nullptr;
    }


    // simulate ===============================================================================================
    Integration::VerletIntegrator verletIntegratorSerial(
            particlePtrsSerial,
            accelerationFunctionSerial, nullptr, nullptr, nullptr, collisionModel);

    Integration::ParallelVerletIntegrator verletIntegratorParallel(
            particlePtrsParallelNew,
            accelerationFunctionParallelNew, nullptr, nullptr, nullptr, collisionModel);

    runIntegrator(verletIntegratorSerial, timeSteps, dt, "serial");
    runIntegrator(verletIntegratorParallel, timeSteps, dt, "parallel");

    std::vector<double> diffMags;
    for (unsigned int i=0; i<nIonsTotal; ++i){
        diffMags.push_back( (particlesSerial[i]->getLocation() - particlesParallelNew[i]->getLocation()).magnitude() );
    }
    double sum = std::accumulate(diffMags.begin(), diffMags.end(), 0.0);

    if (verbose) {
        for (unsigned int i = 0; i<nIonsTotal; ++i) {
            std::cout << particlesSerial[i]->getLocation() << " | " << particlesParallelNew[i]->getLocation() << " | "
                      << (particlesSerial[i]->getLocation()-particlesParallelNew[i]->getLocation()).magnitude()
                      << std::endl;
        }
    }
    std::cout << "sum diff: " << sum <<std::endl;
}