/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2020 - Physical and Theoretical Chemistry /
 Institute of Pure and Applied Mass Spectrometry
 of the University of Wuppertal, Germany

 IDSimF is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 IDSimF is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with IDSimF.  If not, see <https://www.gnu.org/licenses/>.

 ------------
 BT-RS-IMSSim.cpp

 Isothermic continuous field ion mobility spectrometry (IMS) transport and chemistry simulation

 ****************************/


#include "json.h"
#include "parameterParsing.hpp"
#include "RS_Simulation.hpp"
#include "RS_SimulationConfiguration.hpp"
#include "RS_ConfigFileParser.hpp"
#include "RS_ConcentrationFileWriter.hpp"
#include "PSim_verletIntegrator.hpp"
#include "PSim_velocityIntegrator.hpp"
#include "PSim_trajectoryHDF5Writer.hpp"
#include "PSim_util.hpp"
#include "CollisionModel_AbstractCollisionModel.hpp"
#include "CollisionModel_EmptyCollisionModel.hpp"
#include "CollisionModel_HardSphere.hpp"
#include "CollisionModel_StatisticalDiffusion.hpp"
#include "CollisionModel_util.hpp"
#include <iostream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <numeric>
#include <CollisionModel_MultiCollisionModel.hpp>

enum IntegratorType{
    VERLET, SIMPLE, NO_INTEGRATOR
};
enum CollisionModelType{
    HS, SDS, NO_COLLISONS
};

std::string key_ChemicalIndex = "keyChemicalIndex";


int main(int argc, const char *argv[]){

    // open configuration, parse configuration file =========================================
    if (argc < 2) {
        std::cout << "no configuration file given" << std::endl;
        return (0);
    }

    std::string confFileName = argv[1];
    Json::Value confRoot = readConfigurationJson(confFileName);

    std::vector<int> nParticles = intVectorConfParameter("n_particles", confRoot);
    int nSteps = intConfParameter("sim_time_steps", confRoot);
    int concentrationWriteInterval = intConfParameter("concentrations_write_interval", confRoot);
    int trajectoryWriteInterval = intConfParameter("trajectory_write_interval", confRoot);
    bool writeVelocities = boolConfParameter("trajectory_write_velocities", confRoot);
    double dt_s = doubleConfParameter("dt_s", confRoot);
    double eFieldMagnitude = doubleConfParameter("electric_field_mag_Vm-1", confRoot);
    double spaceChargeFactor = doubleConfParameter("space_charge_factor", confRoot);

    double startWidthX_m = doubleConfParameter("start_width_x_mm", confRoot) / 1000.0;
    double startWidthYZ_m = doubleConfParameter("start_width_yz_mm", confRoot) / 1000.0;
    double stopPosX_m = doubleConfParameter("stop_position_x_mm", confRoot) / 1000.0;

    //read and check gas parameters:
    std::string transportModelType = stringConfParameter("transport_model_type", confRoot);
    double backgroundTemperature_K = doubleConfParameter("background_temperature_K", confRoot);

    std::vector<double> backgroundPartialPressures_Pa = doubleVectorConfParameter("background_partial_pressures_Pa", confRoot);
    std::vector<double> collisionGasMasses_Amu = doubleVectorConfParameter("collision_gas_masses_amu", confRoot);
    std::vector<double> collisionGasDiameters_angstrom = doubleVectorConfParameter("collision_gas_diameters_angstrom", confRoot);

    int nBackgroundGases = backgroundPartialPressures_Pa.size();
    if (collisionGasMasses_Amu.size() != nBackgroundGases ||
        collisionGasMasses_Amu.size() != nBackgroundGases ||
        collisionGasDiameters_angstrom.size() != nBackgroundGases)
    {
        throw std::invalid_argument("Inconsistent background gas configuration");
    }

    //compute additional gas parameters:
    double totalBackgroundPressure_Pa = std::accumulate(
            backgroundPartialPressures_Pa.begin(),
            backgroundPartialPressures_Pa.end(), 0.0);

    std::vector<double> collisionGasDiameters_m;
    std::transform(
            collisionGasDiameters_angstrom.begin(),
            collisionGasDiameters_angstrom.end(),
            std::back_inserter(collisionGasDiameters_m),
            [](double cgd)-> double {return cgd*1e-10;} );

    std::string projectFilename;
    if (argc==3) {
        projectFilename = argv[2];
    }
    else {
        std::stringstream ss;
        ss << argv[1] << "_result.txt";
        projectFilename = ss.str();
    }
    // ======================================================================================

    //read and prepare chemical configuration ===============================================
    std::string rsConfFileName = pathRelativeToConfFile(
                                    confFileName,
                                    stringConfParameter("reaction_configuration", confRoot));
    RS::ConfigFileParser parser = RS::ConfigFileParser();
    RS::Simulation rsSim = RS::Simulation(parser.parseFile(rsConfFileName));
    RS::SimulationConfiguration *simConf = rsSim.simulationConfiguration();
    //prepare a map for retrieval of the substance index:
    std::map<RS::Substance *, int> substanceIndices;
    std::vector<RS::Substance *> discreteSubstances = simConf->getAllDiscreteSubstances();
    std::vector<double> ionMobility; // = doubleVectorConfParameter("ion_mobility",confRoot);
    for (int i = 0; i < discreteSubstances.size(); i++) {
        substanceIndices.insert(std::pair<RS::Substance *, int>(discreteSubstances[i], i));
        ionMobility.push_back(discreteSubstances[i]->mobility());
    }

    // prepare file writer  =================================================================
    RS::ConcentrationFileWriter resultFilewriter(projectFilename);

    //prepare auxiliary parameters transform functions
    ParticleSimulation::additionalPartParamFctType additionalParamTFct;
    std::vector<std::string> auxParamNames;

    if (writeVelocities) {
        additionalParamTFct = [](BTree::Particle *particle) -> std::vector<double>
            {
                std::vector<double> result = {
                        particle->getAuxScalarParam(key_ChemicalIndex),
                        particle->getVelocity().x(),
                        particle->getVelocity().y(),
                        particle->getVelocity().z()
                };
                return result;
            };
        auxParamNames = {"chemical_id","velocity x","velocity y","velocity z"};
    }
    else {
        additionalParamTFct = [](BTree::Particle *particle) -> std::vector<double>
        {
            std::vector<double> result = {
                    particle->getAuxScalarParam(key_ChemicalIndex)
            };
            return result;
        };
        auxParamNames = {"chemical_id"};
    }

    //init hdf5 filewriter
    std::string hdf5Filename = projectFilename + "_trajectories.hd5";
    ParticleSimulation::TrajectoryHDF5Writer hdf5Writer(hdf5Filename, auxParamNames, additionalParamTFct);

    int ionsInactive = 0;

    //fixme: nAllParticles and nTotalParticles are the same parameter
    int nAllParticles = 0;
    for (const auto ni: nParticles) {
        nAllParticles += ni;
    }

    // init simulation  =====================================================================

    // create and add simulation particles:
    int nParticlesTotal = 0;
    std::vector<uniqueReactivePartPtr> particles;
    std::vector<BTree::Particle *> particlesPtrs;
    std::vector<std::vector<double>> trajectoryAdditionalParams;

    Core::Vector initCorner(0, 0, 0);
    Core::Vector initBoxSize(startWidthX_m, startWidthYZ_m, startWidthYZ_m);

    for (int i = 0; i < nParticles.size(); i++) {
        RS::Substance *subst = simConf->substance(i);
        std::vector<Core::Vector> initialPositions =
                ParticleSimulation::util::getRandomPositionsInBox(nParticles[i], initCorner, initBoxSize);
        for (int k = 0; k < nParticles[i]; k++) {
            uniqueReactivePartPtr particle = std::make_unique<RS::ReactiveParticle>(subst);

            particle->setLocation(initialPositions[k]);
            particlesPtrs.push_back(particle.get());
            rsSim.addParticle(particle.get(), nParticlesTotal);
            particles.push_back(std::move(particle));
            trajectoryAdditionalParams.push_back(std::vector<double>(1));
            nParticlesTotal++;
        }
    }

    RS::ReactionConditions reactionConditions = RS::ReactionConditions();
    reactionConditions.temperature = backgroundTemperature_K;
    reactionConditions.pressure = totalBackgroundPressure_Pa;
    reactionConditions.electricField = eFieldMagnitude;
    reactionConditions.totalReactionEnergy = 0.0;

    resultFilewriter.initFile(simConf);
    // ======================================================================================

    //check which integrator type we have to setup:
    std::vector<std::string> verletTypes{"btree_SDS", "btree_HS"};
    auto vType = std::find(std::begin(verletTypes), std::end(verletTypes), transportModelType);

    IntegratorType integratorType;
    if (vType != std::end(verletTypes)) {
        integratorType = VERLET;
        std::cout << "Verlet type simulation" << '\n';
    } else if (transportModelType == "simple") {
        integratorType = SIMPLE;
        std::cout << "Simple transport simulation" << '\n';
    } else if (transportModelType == "no_transport") {
        integratorType = NO_INTEGRATOR;
        std::cout << "No transport simulation" << '\n';
    } else {
        throw std::invalid_argument("illegal transport simulation type");
    }
    //===================================================================================================

    // define trajectory integration parameters / functions =================================
    double referencePressure_Pa = 100000;
    double referenceTemperature_K = 273.15;
    double backgroundPTRatio = referencePressure_Pa / totalBackgroundPressure_Pa * backgroundTemperature_K / referenceTemperature_K;

    auto accelerationFctVerlet =
        [eFieldMagnitude, spaceChargeFactor]
        (BTree::Particle *particle, int particleIndex, BTree::Tree &tree, double time, int timestep)
        {
            double particleCharge = particle->getCharge();

            Core::Vector fieldForce(eFieldMagnitude*particleCharge, 0, 0);

            if (spaceChargeFactor==0.0) {
                return (fieldForce/particle->getMass());
            }
            else {
                Core::Vector spaceChargeForce =
                        tree.computeEFieldFromTree(*particle)*(particleCharge*spaceChargeFactor);
                return ((fieldForce+spaceChargeForce)/particle->getMass());
            }
        };


    auto timestepWriteFctSimple =
        [&hdf5Writer, trajectoryWriteInterval, &ionsInactive, nSteps]
        (std::vector<BTree::Particle*> &particles, double time, int timestep, bool lastTimestep)
        {
            if (lastTimestep) {
                hdf5Writer.writeTimestep(particles, time);
                hdf5Writer.writeSplatTimes(particles);
                hdf5Writer.finalizeTrajectory();
                std::cout << "finished ts:" << timestep << " time:" << time << std::endl;
            }
            else if (timestep%trajectoryWriteInterval==0) {
                hdf5Writer.writeTimestep(particles, time);
                std::cout << "ts:" << timestep << " time:" << time << " splatted ions:" << ionsInactive << std::endl;
            }
        };

    auto timestepWriteFctVerlet =
        [&timestepWriteFctSimple]
        (std::vector<BTree::Particle*>& particles, BTree::Tree& tree, double time, int timestep, bool lastTimestep)
        {
            timestepWriteFctSimple(particles, time, timestep, lastTimestep);
        };

    auto otherActionsFunctionIMSSimple =
        [stopPosX_m, &ionsInactive]
        (Core::Vector &newPartPos, BTree::Particle *particle, int particleIndex, double time, int timestep)
        {
            if (newPartPos.x() >= stopPosX_m) {
                particle->setActive(false);
                particle->setSplatTime(time);
                ionsInactive++;
            }
        };

    auto otherActionsFunctionIMSVerlet =
        [&otherActionsFunctionIMSSimple]
        (Core::Vector &newPartPos, BTree::Particle *particle, int particleIndex,
                BTree::Tree &tree, double time, int timestep)
        {
            otherActionsFunctionIMSSimple(newPartPos, particle, particleIndex, time, timestep);
        };


    //define and init transport models and trajectory integrators:
    std::unique_ptr<CollisionModel::AbstractCollisionModel> collisionModelPtr;
    CollisionModelType collisionModelType = NO_COLLISONS;
    if (transportModelType=="btree_SDS") {
        //prepare SDS model
        if (nBackgroundGases!=1) {
            throw std::invalid_argument("SDS simulation requires a single collision gas");
        }

        //create sds collision model, if statistics file is given: create with custom statistics
        std::unique_ptr<CollisionModel::StatisticalDiffusionModel> collisionModel;
        if (isConfFileKey("sds_collision_statistics", confRoot)){
            std::string sdsCollisionStatisticsFileName = pathRelativeToConfFile(
                    confFileName,
                    stringConfParameter("sds_collision_statistics",confRoot));
            std::cout << "SDS with custom collision statistics file: " << sdsCollisionStatisticsFileName << std::endl;
            CollisionModel::CollisionStatistics cs(sdsCollisionStatisticsFileName);
            collisionModel =
                    std::make_unique<CollisionModel::StatisticalDiffusionModel>(
                            backgroundPartialPressures_Pa[0],
                            backgroundTemperature_K,
                            collisionGasMasses_Amu[0],
                            collisionGasDiameters_m[0],
                            cs
                            );
        }
        else {
            collisionModel =
                    std::make_unique<CollisionModel::StatisticalDiffusionModel>(
                            backgroundPartialPressures_Pa[0],
                            backgroundTemperature_K,
                            collisionGasMasses_Amu[0],
                            collisionGasDiameters_m[0]);
        }

        for (const auto& particle: particlesPtrs) {
            collisionModel->setSTPParameters(*particle);
        }
        collisionModelPtr = std::move(collisionModel);
        collisionModelType = SDS;
    }
    else if (transportModelType=="btree_HS") {
        //prepare multimodel with multiple Hard Sphere models (one per collision gas)
        std::vector<std::unique_ptr<CollisionModel::AbstractCollisionModel>> hsModels;
        for (int i = 0; i<nBackgroundGases; ++i) {
            auto hsModel = std::make_unique<CollisionModel::HardSphereModel>(
                    backgroundPartialPressures_Pa[i],
                    backgroundTemperature_K,
                    collisionGasMasses_Amu[i],
                    collisionGasDiameters_m[i]);
            hsModels.emplace_back(std::move(hsModel));
        }

        std::unique_ptr<CollisionModel::MultiCollisionModel> collisionModel =
                std::make_unique<CollisionModel::MultiCollisionModel>(std::move(hsModels));

        collisionModelPtr = std::move(collisionModel);
        collisionModelType = HS;
    }

    //init trajectory simulation object:
    std::unique_ptr<ParticleSimulation::AbstractTimeIntegrator> trajectoryIntegrator;
    if (integratorType == VERLET) {
        trajectoryIntegrator = std::make_unique<ParticleSimulation::VerletIntegrator>(
                particlesPtrs,
                accelerationFctVerlet, timestepWriteFctVerlet, otherActionsFunctionIMSVerlet,
                *collisionModelPtr
        );
    } else if (integratorType == SIMPLE) {

        auto velocityFctSimple = [eFieldMagnitude,backgroundPTRatio](BTree::Particle *particle, int particleIndex, double time, int timestep){
            double particleMobility = particle->getMobility();

            Core::Vector velocity(eFieldMagnitude * particleMobility * backgroundPTRatio, 0, 0);
            return velocity;
        };

        trajectoryIntegrator = std::make_unique<ParticleSimulation::VelocityIntegrator>(
                particlesPtrs,
                velocityFctSimple, timestepWriteFctSimple, otherActionsFunctionIMSSimple
        );
    }
// ======================================================================================


    // simulate   ===========================================================================
    clock_t begin = std::clock();

    for (int step = 0; step < nSteps; step++) {
        if (step % concentrationWriteInterval == 0) {
            rsSim.printConcentrations();
            resultFilewriter.writeTimestep(rsSim);
        }
        for (int i = 0; i < nParticlesTotal; i++) {
            bool reacted = rsSim.react(i, reactionConditions, dt_s);
            int substIndex = substanceIndices.at(particles[i]->getSpecies());
            particles[i]->setAuxScalarParam(key_ChemicalIndex, substIndex);

            if (reacted && collisionModelType == SDS) {
                //we had an reaction event: update the collision model parameters for the particle which are not
                //based on location (mostly STP parameters in SDS)
                collisionModelPtr->initializeModelParameters(*particles[i]);
            }
        }
        rsSim.advanceTimestep(dt_s);
        if (trajectoryIntegrator) {
            trajectoryIntegrator->runSingleStep(dt_s);
        }

        if (ionsInactive >= nAllParticles) {
            break;
        }
    }
    resultFilewriter.writeReactionStatistics(rsSim);
    if (trajectoryIntegrator) {
        trajectoryIntegrator->terminateSimulation();
    }

    clock_t end = std::clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "----------------------"<<std::endl;
    std::cout << "Reaction Events:"<<std::endl;
    rsSim.printReactionStatistics();
    std::cout << "----------------------"<<std::endl;
    std::cout << "total reaction events:" << rsSim.totalReactionEvents() << " ill events:" << rsSim.illEvents() << std::endl;
    std::cout << "ill fraction: " << rsSim.illEvents() / (double) rsSim.totalReactionEvents() << std::endl;

    std::cout << "elapsed seconds " << elapsed_secs << std::endl;

    // ======================================================================================

    return 0;
}
