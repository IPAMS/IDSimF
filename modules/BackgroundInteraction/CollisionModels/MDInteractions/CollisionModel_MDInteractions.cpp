/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2022 - Physical and Theoretical Chemistry /
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
 ****************************/

#include "CollisionModel_MDInteractions.hpp"
#include "Core_math.hpp"
#include "Core_utils.hpp"
#include "Core_randomGenerators.hpp"
#include <cmath>
#include <array>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <functional>
#include <initializer_list>

/**
 * Constructor for static pressure and temperatur.
 * 
 * @param collisionGasPolarizabilityM3 Polarizability of the background gas in m³
 * @param collisionMolecule String identifier of the background gas particle (needs to be contained in the structure map)
 * @param integrationTime Maximum integration time for each collision 
 * @param subTimeStep Timstep length for the trajectory intgeration (leapfrog) or length of first timestep in RKF45
 * @param collisionRadiusScaling Scaling parameter for the collision cross section used for collision probability estimation and defintion of an "actual" collision 
 * @param angleThetaScaling Scaling parameter for the maximum angle under which a collision is still taken as "hit", scaling independently of the collision probability
 * @param spawnRadius Radius of the spawn sphere for the background gas particle 
 * @param molecularStructureCollection Key-value map for all read-in molecular structures which can be used to construct ions and background gas particles
*/
CollisionModel::MDInteractionsModel::MDInteractionsModel(double staticPressure,
                                                        double staticTemperature,
                                                        double collisionGasMassAmu,
                                                        double collisionGasDiameterM,
                                                        std::string collisionMolecule,
                                                        double integrationTime,
                                                        double subTimeStep,
                                                        double collisionRadiusScaling,
                                                        double angleThetaScaling,
                                                        double spawnRadius,
                                                        bool rotationActive,
                                                        std::unique_ptr<CollisionModel::AbstractMDForceField> forceField,
                                                        std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection) :
        MDInteractionsModel(
        getConstantScalarFunction(staticPressure),
        getConstantVectorFunction(Core::Vector(0.0, 0.0, 0.0)),
        staticTemperature,
        collisionGasMassAmu,
        collisionGasDiameterM,
        collisionMolecule,
        integrationTime,
        subTimeStep,
        collisionRadiusScaling,
        angleThetaScaling,
        spawnRadius,
        rotationActive,
        std::move(forceField),
        molecularStructureCollection) { }

CollisionModel::MDInteractionsModel::MDInteractionsModel(double staticPressure,
                                                        double staticTemperature,
                                                        double collisionGasMassAmu,
                                                        double collisionGasDiameterM,
                                                        std::string collisionMolecule,
                                                        double integrationTime,
                                                        double subTimeStep,
                                                        double collisionRadiusScaling,
                                                        double angleThetaScaling,
                                                        double spawnRadius,
                                                        std::unique_ptr<CollisionModel::AbstractMDForceField> forceField,
                                                        std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection) :
        MDInteractionsModel(
        getConstantScalarFunction(staticPressure),
        getConstantVectorFunction(Core::Vector(0.0, 0.0, 0.0)),
        staticTemperature,
        collisionGasMassAmu,
        collisionGasDiameterM,
        collisionMolecule,
        integrationTime,
        subTimeStep,
        collisionRadiusScaling,
        angleThetaScaling,
        spawnRadius,
        false,
        std::move(forceField),
        molecularStructureCollection) { }

CollisionModel::MDInteractionsModel::MDInteractionsModel(std::function<double(Core::Vector& location)> pressureFunction,
                                                        std::function<Core::Vector(Core::Vector& location)> velocityFunction,
                                                        double staticTemperature,
                                                        double collisionGasMassAmu,
                                                        double collisionGasDiameterM,
                                                        std::string collisionMolecule,
                                                        double integrationTime,
                                                        double subTimeStep,
                                                        double collisionRadiusScaling,
                                                        double angleThetaScaling,
                                                        double spawnRadius,
                                                        bool rotationActive,
                                                        std::unique_ptr<CollisionModel::AbstractMDForceField> forceField,
                                                        std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection) :
        MDInteractionsModel(
                std::move(pressureFunction),
                std::move(velocityFunction),
                getConstantScalarFunction(staticTemperature),
                collisionGasMassAmu,
                collisionGasDiameterM,
                collisionMolecule,
                integrationTime,
                subTimeStep,
                collisionRadiusScaling,
                angleThetaScaling,
                spawnRadius,
                rotationActive,
                std::move(forceField),
                molecularStructureCollection) { }

CollisionModel::MDInteractionsModel::MDInteractionsModel(std::function<double(Core::Vector& location)> pressureFunction,
                                                        std::function<Core::Vector(Core::Vector& location)> velocityFunction,
                                                        std::function<double(const Core::Vector&)> temperatureFunction,
                                                        double collisionGasMassAmu,
                                                        double collisionGasDiameterM,
                                                        std::string collisionMolecule,
                                                        double integrationTime,
                                                        double subTimeStep,
                                                        double collisionRadiusScaling,
                                                        double angleThetaScaling,
                                                        double spawnRadius,
                                                        std::unique_ptr<CollisionModel::AbstractMDForceField> forceField,
                                                        std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection) :
        MDInteractionsModel(
                std::move(pressureFunction),
                std::move(velocityFunction),
                std::move(temperatureFunction),
                collisionGasMassAmu,
                collisionGasDiameterM,
                collisionMolecule,
                integrationTime,
                subTimeStep,
                collisionRadiusScaling,
                angleThetaScaling,
                spawnRadius,
                false,
                std::move(forceField),
                molecularStructureCollection) { }

CollisionModel::MDInteractionsModel::MDInteractionsModel(std::function<double(Core::Vector& location)> pressureFunction,
                                                        std::function<Core::Vector(Core::Vector& location)> velocityFunction,
                                                        std::function<double(const Core::Vector&)> temperatureFunction,
                                                        double collisionGasMassAmu,
                                                        double collisionGasDiameterM,
                                                        std::string collisionMolecule,
                                                        double integrationTime,
                                                        double subTimeStep,
                                                        double collisionRadiusScaling,
                                                        double angleThetaScaling,
                                                        double spawnRadius,
                                                        bool rotationActive,
                                                        std::unique_ptr<CollisionModel::AbstractMDForceField> forceField,
                                                        std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection) :
        MDIntegrator(collisionMolecule, rotationActive, std::move(forceField)),
        pressureFunction_(std::move(pressureFunction)),
        velocityFunction_(std::move(velocityFunction)),
        temperatureFunction_(std::move(temperatureFunction)),
        collisionGasMass_kg_(collisionGasMassAmu*Core::AMU_TO_KG),
        collisionGasDiameter_m_(collisionGasDiameterM),
        integrationTime_(integrationTime),
        subTimeStep_(subTimeStep),
        collisionRadiusScaling_(collisionRadiusScaling),
        angleThetaScaling_(angleThetaScaling),
        spawnRadius_(spawnRadius),
        molecularStructureCollection_(std::move(molecularStructureCollection))
{}



void CollisionModel::MDInteractionsModel::initializeModelParticleParameters(Core::Particle& /*ion*/) const {

}

void CollisionModel::MDInteractionsModel::updateModelParticleParameters(Core::Particle& /*ion*/) const {

}

/**
 * Updates trajectory recording if timestep recording parameter is exceeded
*/
void CollisionModel::MDInteractionsModel::updateModelTimestepParameters(unsigned int timestep, double /*time*/) {
    if (legacyTWriterConf_.modelRecordsTrajectory && timestep >= recordTrajectoryStartTimeStep_){
        legacyTWriterConf_.recordingActive = true;
    }
    if (hdf5TWriterConf_.modelRecordsTrajectory && timestep >= recordTrajectoryStartTimeStep_){
        hdf5TWriterConf_.recordingActive = true;
    }
}

void CollisionModel::MDInteractionsModel::modifyAcceleration(Core::Vector& /*acceleration*/, Core::Particle& /*particle*/,
                                                         double /*dt*/) {

}

/**
 * Modifies the velocity of the particle based on a molecular dynamics approach.
 * Collision probability is estimated by a hard-sphere model.
 * Trajectory is checked for energy conservation of 10 % and if necessary is repeated for up to 
 * 100 times under modification of the starting timestep length. 
 * @param particle particle whose velocity is to be modified 
 * @param dt timestep length of overarching ion simulation 
 */
void CollisionModel::MDInteractionsModel::modifyVelocity(Core::Particle& particle, double dt) {
    Core::RandomSource* rndSource = Core::globalRandomGeneratorPool->getThreadRandomSource();

    // Calculate collision cross section between particle and collision gas:
    double collisionRadius = collisionRadiusScaling_*(particle.getDiameter() + collisionGasDiameter_m_)/2.0;
    double sigma_m2 = M_PI * collisionRadius * collisionRadius;

    Core::Vector moleculeComPosition = particle.getLocation();
    double localPressure_Pa = pressureFunction_(moleculeComPosition);
    if (Core::isDoubleEqual(localPressure_Pa, 0.0)){
        return; //pressure 0 means no collision at all
    }

    // Transform the frame of reference in a frame where the mean background gas velocity is zero.
    Core::Vector vGasMean = velocityFunction_(moleculeComPosition);
    Core::Vector vFrameMeanBackRest = particle.getVelocity() - vGasMean;

    double vRelIonMeanBackRest = vFrameMeanBackRest.magnitude(); //relative ion relative to bulk gas velocity

    // Calculate the mean free path (MFP) from current ion velocity:

    // a static ion leads in static gas leads to a relative velocity of zero, which leads
    // to undefined behavior due to division by zero later.
    // The whole process converges to the MFP and collision probability of a static ion, thus
    // it is possible to assume a small velocity (1 nm/s) for the static ions to get rid of undefined behavior
    if (vRelIonMeanBackRest < 1e-9){
        vRelIonMeanBackRest = 1e-9;
    }

    // Calculate the mean gas speed (m/s)
    double temperature_K = temperatureFunction_(moleculeComPosition);
    double vMeanGas = std::sqrt(8.0*Core::K_BOLTZMANN*temperature_K/M_PI/(collisionGasMass_kg_));

    // Calculate the median gas speed (m/s)
    double vMedianGas = std::sqrt(2.0*Core::K_BOLTZMANN*temperature_K/(collisionGasMass_kg_));

    // Compute the mean relative speed (m/s) between ion and gas.
    double s = vRelIonMeanBackRest / vMedianGas;
    double cMeanRel = vMeanGas * (
            (s + 1.0/(2.0*s)) * 0.5 * sqrt(M_PI) * std::erf(s) + 0.5 * std::exp(-s*s) );

    // Compute mean-free-path (m)
    double effectiveMFP_m = Core::K_BOLTZMANN * temperature_K *
                            (vRelIonMeanBackRest / cMeanRel) / (localPressure_Pa * sigma_m2);

    // Compute probability of collision in the current time-step.
    double collisionProb = 1.0 - std::exp(-vRelIonMeanBackRest * dt / effectiveMFP_m);

    // FIXME: The time step length dt is unrestricted
    // Possible mitigation: Throw warning / exception if collision probability becomes too high
    if(collisionProb > 0.20)
        std::cout << "collisionProb " << collisionProb << '\n';
    // Decide if a collision actually happens:
    if (rndSource->uniformRealRndValue() > collisionProb){
        return; // no collision takes place
    }
    
    bool trajectorySuccess = false;
    int iterations = 0;
    double spawnRad = spawnRadius_;
    double collisionTheta = std::asin(collisionRadius / spawnRad);
    double tolerance = 1e-8;


    do{
        // Collision happens
        // Construct the actual molecule and its atoms
        CollisionModel::Molecule mole = CollisionModel::Molecule(Core::Vector(0.0, 0.0, 0.0), Core::Vector(0.0, 0.0, 0.0), particle.getMolecularStructure());

        // Construct the background gas particle
        CollisionModel::Molecule bgMole = CollisionModel::Molecule(Core::Vector(0.0, 0.0, 0.0), Core::Vector(0.0, 0.0, 0.0),
                                            molecularStructureCollection_.at(collisionMolecule_));

        //Init new HDF5 trajectory
        if(hdf5TWriterConf_.recordingActive == true){
            hdf5TrajectoryWriter_->initNewTrajectory(mole, bgMole);
        }


        // Give background gas its position, velocity, rotation:
        // Calculate the standard deviation of the one dimensional velocity distribution of the
        // background gas particles. Std. dev. in one dimension is given from Maxwell-Boltzmann
        // as sqrt(kT / particle mass).
        double  vrStdevBgMolecule = std::sqrt( Core::K_BOLTZMANN * temperature_K / (collisionGasMass_kg_) );
        Core::Vector velocityBgMolecule = { rndSource->normalRealRndValue() * vrStdevBgMolecule - particle.getVelocity().x(),
                                            rndSource->normalRealRndValue() * vrStdevBgMolecule - particle.getVelocity().y(),
                                            rndSource->normalRealRndValue() * vrStdevBgMolecule - particle.getVelocity().z()};

        bgMole.setComVel(velocityBgMolecule);

        // calculate random point on sphere
        // as follows:
        // draw random number in as long until magnitude is less than 1
        // normalize result
        double circleVectorMagnitude = 0;
        Core::Vector circleVector(0,0,0);
        double directionAngle = 0;
        do{
            circleVector = Core::Vector{
                  (rndSource->uniformRealRndValue() * 2 - 1)
                , (rndSource->uniformRealRndValue() * 2 - 1)
                , (rndSource->uniformRealRndValue() * 2 - 1)
            };
            circleVectorMagnitude = circleVector.magnitude();
            circleVector = spawnRad / circleVectorMagnitude * circleVector;
            directionAngle  = std::acos(
                ( (-1. * circleVector) * velocityBgMolecule)
                / ( spawnRad * velocityBgMolecule.magnitude() )
            );
        }while(circleVectorMagnitude > 1 || directionAngle >  angleThetaScaling_ * collisionTheta);

        bgMole.setComPos(circleVector);

        double energyRotMolecule = 0;
        double energyRotBg = 0;
        if(rotationActive_){
            energyRotMolecule = initRotation(mole, temperature_K);
            energyRotBg = initRotation(bgMole, temperature_K);
        }

        std::vector<CollisionModel::Molecule*> moleculesPtr = {&mole, &bgMole};

        for(auto* molecule : moleculesPtr){
            Core::Matrix3 initialRotationMatrix = molecule->calcRotationMatrix(rndSource->uniformRealRndValue()*2*M_PI-M_PI, 
                                                                               rndSource->uniformRealRndValue()*2*M_PI-M_PI, 
                                                                               rndSource->uniformRealRndValue()*2*M_PI-M_PI);
            molecule->setRotationMatrix(initialRotationMatrix);
            //molecule->rotateMoleculeRotationMatrix();
        }

        // possible check for energy conservation
        std::vector<Core::Vector> startVelocity;
        double kineticEnergyStart = 0, rotationEnergyStart = 0;

        double startEnergy = 0;
        for(auto* molecule : moleculesPtr){
            startVelocity.push_back(molecule->getComVel());
            kineticEnergyStart += 0.5 * molecule->getMass() * molecule->getComVel().magnitudeSquared();
        }
        if(rotationActive_){
            rotationEnergyStart += (energyRotMolecule + energyRotBg);
        }
        startEnergy = kineticEnergyStart + rotationEnergyStart;

        // Call the sub-integrator
        double finalTime = integrationTime_; //  final integration time in seconds
        double timeStep = subTimeStep_; // step size in seconds

        trajectorySuccess = rk4InternAdaptiveStep(moleculesPtr, timeStep, finalTime, 1e10, collisionRadius, tolerance);
        //trajectorySuccess = leapfrogIntern(moleculesPtr, timeStep, finalTime, collisionRadius);
        double kineticEnergyEnd = 0, rotationEnergyEnd = 0;
        double endEnergy = 0;
        for(auto* molecule : moleculesPtr){
            kineticEnergyEnd += 0.5 * molecule->getMass() * molecule->getComVel().magnitudeSquared();
            if(rotationActive_) rotationEnergyEnd += calcRotEnergy(molecule->getAngVel(), molecule->getInertiaMatrix());
        }

        // check if energy is conserved up to 10% 
        // if not halve the starting timestep length
        endEnergy = kineticEnergyEnd + rotationEnergyEnd; 
        if(!rotationActive_){
            if(startEnergy >= endEnergy*1.10  || startEnergy <= endEnergy*0.90){
                std::cout << "Energy not conserved: " << startEnergy << " " << endEnergy << std::endl;
                trajectorySuccess = false;
                dt = dt/2;
            }
        }else{
            if(startEnergy >= endEnergy*1.10  || startEnergy <= endEnergy*0.90){
                std::cout << "Energy not conserved: " << startEnergy << " " << endEnergy << std::endl;
                std::cout << "Kinetic Energy: " << kineticEnergyStart << " " << kineticEnergyEnd << std::endl;
                std::cout << "Rotation Energy: " << rotationEnergyStart << " " << rotationEnergyEnd << std::endl;
                trajectorySuccess = false;
            }
        }

        if(trajectorySuccess){
            particle.setVelocity(mole.getComVel() + particle.getVelocity() + vGasMean);
        }
        ++iterations;
    }while(!trajectorySuccess && iterations < 100);

    if(trajectorySuccess == false){
        std::cerr << "No trajectory that hit the collision sphere was found or energy could not be conserved.\n";
    }
}

void CollisionModel::MDInteractionsModel::modifyPosition(Core::Particle& /*particle*/, double /*dt*/) {}

