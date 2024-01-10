/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2024 - Physical and Theoretical Chemistry /
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
 test_MDVerlet.cpp

 ****************************/

#include "Integration_velocityIntegrator.hpp"
#include "Integration_verletIntegrator.hpp"
#include "Core_vector.hpp"
#include "Core_particle.hpp"
#include "catch.hpp"
#include "CollisionModel_MDInteractions.hpp"
#include "CollisionModel_MDForceField_LJ12_6.hpp"
#include "FileIO_MolecularStructureReader.hpp"
#include "FileIO_trajectoryHDF5Writer.hpp"

#include "PSim_util.hpp"
#include "PSim_constants.hpp"
#include "CollisionModel_AbstractCollisionModel.hpp"
#include "CollisionModel_HardSphere.hpp"
#include "CollisionModel_StatisticalDiffusion.hpp"
#include "CollisionModel_MultiCollisionModel.hpp"
#include "appUtils_simulationConfiguration.hpp"
#include <iostream>


std::string key_ChemicalIndex = "keyChemicalIndex";

TEST_CASE( "Test MD and integrator", "[ParticleSimulation][VelocityIntegrator][trajectory integration]") {

    SECTION( "Integrate") {
        
        std::string mdCollisionConfFile = "test_molecularstructure_reader.json";
        FileIO::MolecularStructureReader mdConfReader = FileIO::MolecularStructureReader();
        std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection = mdConfReader.readMolecularStructure(mdCollisionConfFile);

        FileIO::partAttribTransformFctType additionalParamTFct;
        std::vector<std::string> auxParamNames;

        additionalParamTFct = [](Core::Particle* particle) -> std::vector<double> {
            std::vector<double> result = {
                    particle->getVelocity().x()
            };
            return result;
        };
        auxParamNames = {"velocity x"};
        

        //init hdf5 filewriter
        std::string hdf5Filename = "test_trajectories.hd5";
        FileIO::TrajectoryHDF5Writer hdf5Writer(hdf5Filename);
        hdf5Writer.setParticleAttributes(auxParamNames, additionalParamTFct);


        //Test with verlet integration:
        double dt = 4e-10;
        unsigned int nParticles = 1;
        int timeSteps = 1;

        //double ionVelocity = 10.0;

        // auto velocityFct = [ionVelocity](Core::Particle* /*particle*/, int /*particleIndex*/, double /*time*/, unsigned int /*timestep*/){
        //     Core::Vector result(ionVelocity, 0, ionVelocity*0.1);
        //     return (result);
        // };

        unsigned int ionsInactive = 0;
        // double referencePressure_Pa = 100000;
        // double referenceTemperature_K = 273.15;
        
        auto accelerationFctVerlet =
                []
                        (Core::Particle* particle, int /*particleIndex*/, SpaceCharge::FieldCalculator& /*scFieldCalculator*/, double /*time*/, int /*timestep*/) {
                    Core::Vector fieldForce(2e-15, 0, 0);
                    return (fieldForce/particle->getMass());

                };

        auto timestepWriteFctSimple =
                [&hdf5Writer]
                        (std::vector<Core::Particle*>& particles, double time, int timestep, bool lastTimestep) {
                    if (lastTimestep) {
                        hdf5Writer.writeTimestep(particles, time);
                        hdf5Writer.writeSplatTimes(particles);
                        hdf5Writer.finalizeTrajectory();
                    }
                    else if (timestep%1==0) {
                        hdf5Writer.writeTimestep(particles, time);
                    }
                };

        auto timestepWriteFctVerlet =
                [&timestepWriteFctSimple]
                        (std::vector<Core::Particle*>& particles,  double time, int timestep,
                         bool lastTimestep) {
                    timestepWriteFctSimple(particles, time, timestep, lastTimestep);
                };

        
        auto otherActionsFunctionIMSSimple =
                [&ionsInactive]
                        (Core::Vector& newPartPos, Core::Particle* particle, int /*particleIndex*/, double time,
                         int /*timestep*/) {
                    if (newPartPos.x()>=1) {
                        particle->setActive(false);
                        particle->setSplatTime(time);
                        ionsInactive++;
                    }
                };
        
        auto otherActionsFunctionIMSVerlet =
                [&otherActionsFunctionIMSSimple]
                        (Core::Vector& newPartPos, Core::Particle* particle, int particleIndex,
                          double time, int timestep) {
                    otherActionsFunctionIMSSimple(newPartPos, particle, particleIndex, time, timestep);
                };


        std::vector<Core::uniquePartPtr>particles;
        std::vector<Core::Particle*>particlesPtrs;

        double yPos = 0;
        for (unsigned int i=0; i<nParticles; ++i){
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
            yPos = yPos+0.01;
        }

        std::unique_ptr<CollisionModel::AbstractCollisionModel> collisionModelPtr;
        std::vector<std::unique_ptr<CollisionModel::AbstractCollisionModel>> mdModels;
        CollisionModel::MDForceField_LJ12_6 forceField(0.203E-30);
        auto forceFieldPtr = std::make_unique<CollisionModel::MDForceField_LJ12_6>(forceField);
        auto mdModel = std::make_unique<CollisionModel::MDInteractionsModel>(
                        100000,
                        298,
                        4.003,
                        2.89e-10,
                        "He",
                        400e-14, 
                        1e-15, 
                        2,
                        4,
                        25e-10,
                        std::move(forceFieldPtr),
                        molecularStructureCollection);
        mdModels.emplace_back(std::move(mdModel));
        std::unique_ptr<CollisionModel::MultiCollisionModel> collisionModel =
                    std::make_unique<CollisionModel::MultiCollisionModel>(std::move(mdModels));
        collisionModelPtr = std::move(collisionModel);

        std::unique_ptr<Integration::AbstractTimeIntegrator> trajectoryIntegrator = std::make_unique<Integration::VerletIntegrator>(
                    particlesPtrs,
                    accelerationFctVerlet, timestepWriteFctVerlet, otherActionsFunctionIMSVerlet,
                    ParticleSimulation::noFunction,
                    collisionModelPtr.get());
        for (int step = 0; step<timeSteps; step++) {
            trajectoryIntegrator->runSingleStep(dt);
        }

        hdf5Writer.finalizeTrajectory();
    }
}