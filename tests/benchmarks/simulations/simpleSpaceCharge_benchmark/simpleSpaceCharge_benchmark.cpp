#include "BTree_node.hpp"
#include "BTree_parallelNodeOriginal.hpp"
#include "BTree_tree.hpp"
#include "BTree_parallelTree.hpp"
#include "BTree_parallelTreeOriginal.hpp"
#include "BTree_particle.hpp"
#include "PSim_verletIntegrator.hpp"
#include "PSim_parallelVerletIntegrator.hpp"
#include "PSim_trajectoryHDF5Writer.hpp"
#include "PSim_util.hpp"
#include "CollisionModel_StatisticalDiffusion.hpp"
#include "appUtils_stopwatch.hpp"
#include <iostream>
#include <numeric>

void runIntegrator(ParticleSimulation::AbstractTimeIntegrator &integrator, int timeSteps, double dt, std::string message){
    std::cout << "Benchmark " <<message << std::endl;
    AppUtils::Stopwatch stopWatch;
    stopWatch.start();

    integrator.run(timeSteps, dt);

    stopWatch.stop();
    std::cout << "elapsed wall time:"<< stopWatch.elapsedSecondsWall()<<std::endl;
    std::cout << "elapsed cpu time:"<< stopWatch.elapsedSecondsCPU()<<std::endl;

}


int prepareIons(std::vector<std::unique_ptr<BTree::Particle>> &particles,
                 std::vector<BTree::Particle*> &particlePtrs, int nPerDirection){

    int nTotal = 0;
    for (int i=0; i<nPerDirection; i++){
        double pX = i*1.0/nPerDirection;
        for (int j=0; j<nPerDirection; j++){
            double pY = j*1.0/nPerDirection;
            for (int k=0; k<nPerDirection; k++){
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

int main(int argc, const char * argv[]) {

    bool useCollisionModel = false;
    bool verbose = false;
    int nIonsPerDirection = 23;
    int timeSteps = 200;
    int trajectoryWriteInterval = 10;
    double dt = 1e-3;
    double spaceChargeFactor = 1.0;

    // define functions for the trajectory integration ==================================================
    auto accelerationFunctionSerial =
            [spaceChargeFactor](
                    BTree::Particle *particle, int particleIndex,
                    BTree::Tree &tree, double time, int timestep) -> Core::Vector{

                Core::Vector pos = particle->getLocation();
                double particleCharge = particle->getCharge();

                Core::Vector spaceChargeForce(0,0,0);
                if (spaceChargeFactor > 0) {
                    spaceChargeForce =
                            tree.computeEFieldFromTree(*particle) * (particleCharge * spaceChargeFactor);
                }
                return (spaceChargeForce / particle->getMass());
            };


    auto accelerationFunctionParallel =
            [spaceChargeFactor](
                    BTree::Particle *particle, int particleIndex,
                    BTree::ParallelTreeOriginal &tree, double time, int timestep,std::vector<BTree::ParallelNodeOriginal*> &MyNodes) -> Core::Vector{

                Core::Vector pos = particle->getLocation();
                double particleCharge = particle->getCharge();

                Core::Vector spaceChargeForce(0,0,0);
                if (spaceChargeFactor > 0) {
                    spaceChargeForce =
                            tree.computeEFieldFromTree(*particle,MyNodes) * (particleCharge * spaceChargeFactor);
                }
                return (spaceChargeForce / particle->getMass());
            };

    auto accelerationFunctionParallelNew =
            [spaceChargeFactor](
                    BTree::Particle *particle, int particleIndex,
                    BTree::ParallelTree &tree, double time, int timestep) -> Core::Vector{

                double particleCharge = particle->getCharge();

                Core::Vector spaceChargeForce(0,0,0);
                if (spaceChargeFactor > 0) {
                    spaceChargeForce =
                            tree.computeEFieldFromTree(*particle) * (particleCharge * spaceChargeFactor);
                }
                return (spaceChargeForce / particle->getMass());
            };

    auto hdf5Writer = std::make_unique<ParticleSimulation::TrajectoryHDF5Writer>(
            "test_trajectories.hd5");

    auto timestepWriteFunctionSerialWriting =
            [trajectoryWriteInterval, &hdf5Writer]( std::vector<BTree::Particle *> &particles, BTree::Tree &tree, double time,
                int timestep, bool lastTimestep) {

            if (lastTimestep) {
                hdf5Writer->writeTimestep(particles,time);

                hdf5Writer->writeSplatTimes(particles);
                hdf5Writer->finalizeTrajectory();
                std::cout << "finished ts:" << timestep << " time:" << time << std::endl;
            }

            else if (timestep % trajectoryWriteInterval == 0) {

                std::cout << "ts:" << timestep << " time:" << time << std::endl;
                hdf5Writer->writeTimestep(particles,time);
            }
    };


    std::vector<std::unique_ptr<BTree::Particle>> particlesSerial;
    std::vector<BTree::Particle*>particlePtrsSerial;
    std::vector<std::unique_ptr<BTree::Particle>> particlesParallel;
    std::vector<BTree::Particle*>particlePtrsParallel;
    std::vector<std::unique_ptr<BTree::Particle>> particlesParallelNew;
    std::vector<BTree::Particle*>particlePtrsParallelNew;

    int nIonsTotal = prepareIons(particlesSerial, particlePtrsSerial, nIonsPerDirection);
    prepareIons(particlesParallel, particlePtrsParallel, nIonsPerDirection);
    prepareIons(particlesParallelNew, particlePtrsParallelNew, nIonsPerDirection);
    
    CollisionModel::StatisticalDiffusionModel sdsCollisionModel(100000.0, 298, 28, 3.64e-9);
    CollisionModel::AbstractCollisionModel *collisionModel;
    if (useCollisionModel){
        collisionModel = &sdsCollisionModel;
        for (int i=0; i<nIonsTotal; ++i){
            sdsCollisionModel.setSTPParameters(*particlePtrsSerial[i]);
            sdsCollisionModel.setSTPParameters(*particlePtrsParallel[i]);
            sdsCollisionModel.setSTPParameters(*particlePtrsParallelNew[i]);
        }

    }
    else {
        collisionModel = nullptr;
    }


    // simulate ===============================================================================================
    ParticleSimulation::VerletIntegrator verletIntegratorSerial(
            particlePtrsSerial,
            accelerationFunctionSerial, nullptr, nullptr, nullptr, collisionModel);

    ParticleSimulation::ParallelVerletIntegrator verletIntegratorParallel(
            particlePtrsParallelNew,
            accelerationFunctionParallelNew, nullptr, nullptr, nullptr, collisionModel);

    runIntegrator(verletIntegratorSerial, timeSteps, dt, "serial");
    runIntegrator(verletIntegratorParallel, timeSteps, dt, "parallel");

    std::vector<double> diffMags;
    std::vector<double> diffMagsParallel;
    for (int i=0; i<nIonsTotal; ++i){
        diffMags.push_back( (particlesSerial[i]->getLocation() - particlesParallel[i]->getLocation()).magnitude() );
        diffMagsParallel.push_back( (particlesParallel[i]->getLocation() - particlesParallelNew[i]->getLocation()).magnitude() );
    }
    double sum = std::accumulate(diffMags.begin(), diffMags.end(), 0.0);
    double maximumDiff = *std::max_element(diffMags.begin(), diffMags.end());

    double sumParallel = std::accumulate(diffMagsParallel.begin(), diffMagsParallel.end(), 0.0);
    double maximumDiffParallel = *std::max_element(diffMagsParallel.begin(), diffMagsParallel.end());

    if (verbose) {
        for (int i = 0; i<nIonsTotal; ++i) {
            std::cout << particlesParallel[i]->getLocation() << " | " << particlesParallelNew[i]->getLocation() << " | "
                      << (particlesParallel[i]->getLocation()-particlesParallelNew[i]->getLocation()).magnitude()
                      << std::endl;
        }
    }

    std::cout << "sum diff: " << sum <<std::endl;
    std::cout << "sum diff parallel: " << sumParallel <<std::endl;
}