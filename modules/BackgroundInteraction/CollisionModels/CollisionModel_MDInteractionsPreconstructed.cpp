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

#include "CollisionModel_MDInteractionsPreconstructed.hpp"
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
#include <limits>

CollisionModel::MDInteractionsModelPreconstructed::MDInteractionsModelPreconstructed(double staticPressure,
                                                        double staticTemperature,
                                                        double collisionGasMassAmu,
                                                        double collisionGasDiameterM, 
                                                        double collisionGasPolarizabilityM3,
                                                        std::string collisionMolecule,
                                                        double integrationTime,
                                                        double subTimeStep,
                                                        double collisionRadiusScaling,
                                                        double angleThetaScaling,
                                                        double spawnRadius,
                                                        std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection) :
        MDInteractionsModelPreconstructed(
        getConstantDoubleFunction(staticPressure),
        getConstantVectorFunction(Core::Vector(0.0, 0.0, 0.0)),
        staticTemperature,
        collisionGasMassAmu,
        collisionGasDiameterM,
        collisionGasPolarizabilityM3,
        collisionMolecule,
        integrationTime,
        subTimeStep,
        collisionRadiusScaling,
        angleThetaScaling,
        spawnRadius,
        molecularStructureCollection) { }

CollisionModel::MDInteractionsModelPreconstructed::MDInteractionsModelPreconstructed(std::function<double(Core::Vector& location)> pressureFunction,
                                                        std::function<Core::Vector(Core::Vector& location)> velocityFunction,
                                                        double staticTemperature,
                                                        double collisionGasMassAmu,
                                                        double collisionGasDiameterM,
                                                        double collisionGasPolarizabilityM3,
                                                        std::string collisionMolecule,
                                                        double integrationTime,
                                                        double subTimeStep,
                                                        double collisionRadiusScaling,
                                                        double angleThetaScaling,
                                                        double spawnRadius,
                                                        std::unordered_map<std::string,
                                                        std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection) :
        MDInteractionsModelPreconstructed(
                std::move(pressureFunction),
                std::move(velocityFunction),
                getConstantDoubleFunction(staticTemperature),
                collisionGasMassAmu,
                collisionGasDiameterM,
                collisionGasPolarizabilityM3,
                collisionMolecule,
                integrationTime,
                subTimeStep,
                collisionRadiusScaling,
                angleThetaScaling,
                spawnRadius,
                molecularStructureCollection) { }

CollisionModel::MDInteractionsModelPreconstructed::MDInteractionsModelPreconstructed(std::function<double(Core::Vector& location)> pressureFunction,
                                                        std::function<Core::Vector(Core::Vector& location)> velocityFunction,
                                                        std::function<double(const Core::Vector&)> temperatureFunction,
                                                        double collisionGasMassAmu,
                                                        double collisionGasDiameterM,
                                                        double collisionGasPolarizabilityM3,
                                                        std::string collisionMolecule,
                                                        double integrationTime,
                                                        double subTimeStep,
                                                        double collisionRadiusScaling,
                                                        double angleThetaScaling,
                                                        double spawnRadius,
                                                        std::unordered_map<std::string,
                                                        std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection) :

        collisionGasMass_kg_(collisionGasMassAmu*Core::AMU_TO_KG),
        collisionGasDiameter_m_(collisionGasDiameterM),
        collisionGasPolarizability_m3_(collisionGasPolarizabilityM3),
        collisionMolecule_(collisionMolecule),
        integrationTime_(integrationTime),
        subTimeStep_(subTimeStep),
        collisionRadiusScaling_(collisionRadiusScaling),
        angleThetaScaling_(angleThetaScaling),
        spawnRadius_(spawnRadius),
        pressureFunction_(std::move(pressureFunction)),
        velocityFunction_(std::move(velocityFunction)),
        temperatureFunction_(std::move(temperatureFunction)),
        molecularStructureCollection_(std::move(molecularStructureCollection)) {}

/**
 * Activates trajectory writing and sets trajectory writer configuration
 * @param trajectoryFileName
 * @param trajectoryDistance
 */
void CollisionModel::MDInteractionsModelPreconstructed::setTrajectoryWriter(const std::string& trajectoryFileName,
                                                              double trajectoryDistance,
                                                              int recordTrajectoryStartTimestep) {

    trajectoryOutputStream_ = std::make_unique<std::ofstream>();
    trajectoryOutputStream_->open(trajectoryFileName);

    if (trajectoryOutputStream_->good()){
        recordTrajectoryStartTimeStep_ = recordTrajectoryStartTimestep;
        trajectoryDistance_ = trajectoryDistance;
        modelRecordsTrajectories_ = true;
    }
    else{
        throw (std::runtime_error("Trajectory Output Stream failed to open"));
    }
}

double CollisionModel::MDInteractionsModelPreconstructed::calcSign(double value){
    if(value > 0){
        return 1.;
    }else if(value < 0){
        return -1.;
    }else{
        return 0;
    }
}

void CollisionModel::MDInteractionsModelPreconstructed::writeTrajectory(double distance, Core::Vector positionBgMolecule, Core::Vector velocityBgMolecule, 
                        std::vector<Core::Vector> forceMolecules, bool endOfTrajectory, std::ofstream* file, double time, double dt){
    if(distance < trajectoryDistance_){
        *file << positionBgMolecule.x() << ", " << positionBgMolecule.y() << ", " << positionBgMolecule.z() << 
        ", " << distance << ", " << time <<
        ", " << velocityBgMolecule.x() << ", " << velocityBgMolecule.y() << ", " << velocityBgMolecule.z() << 
        ", " << forceMolecules[1].x() << ", " << forceMolecules[1].y() << ", " << forceMolecules[1].z() <<
        ", " << dt << 
        std::endl;
    }
    if(endOfTrajectory == true){
        *file << "###" << std::endl;
    }
}

void CollisionModel::MDInteractionsModelPreconstructed::initializeModelParticleParameters(Core::Particle& /*ion*/) const {

}

void CollisionModel::MDInteractionsModelPreconstructed::updateModelParticleParameters(Core::Particle& /*ion*/) const {

}

void CollisionModel::MDInteractionsModelPreconstructed::updateModelTimestepParameters(int timestep, double /*time*/) {
    if (modelRecordsTrajectories_ && timestep > recordTrajectoryStartTimeStep_){
        trajectoryRecordingActive_ = true;
    }
}

void CollisionModel::MDInteractionsModelPreconstructed::modifyAcceleration(Core::Vector& /*acceleration*/, Core::Particle& /*particle*/,
                                                         double /*dt*/) {

}

void CollisionModel::MDInteractionsModelPreconstructed::modifyVelocity(Core::Particle& particle, double dt) {
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

    // Always let collision happen  

    bool trajectorySuccess = false;
    int iterations = 0;
    double spawnRad = spawnRadius_;
    double collisionTheta = std::asin(collisionRadius / spawnRad);


    do{
        // Collision happens
        // Construct the actual molecule and its atoms
        CollisionModel::Molecule mole = CollisionModel::Molecule(Core::Vector(0.0, 0.0, 0.0), Core::Vector(0.0, 0.0, 0.0), particle.getMolecularStructure());

        // Construct the background gas particle
        CollisionModel::Molecule bgMole = CollisionModel::Molecule(Core::Vector(0.0, 0.0, 0.0), Core::Vector(0.0, 0.0, 0.0),
                                            molecularStructureCollection_.at(collisionMolecule_));

        // Give background gas its position, velocity, rotation:
        // Calculate the standard deviation of the one dimensional velocity distribution of the
        // background gas particles. Std. dev. in one dimension is given from Maxwell-Boltzmann
        // as sqrt(kT / particle mass).

        // HEAD-ON COLLISION IN X 
        bgMole.setComPos(Core::Vector(spawnRad, 0, 0));
        double  vrStdevBgMolecule = std::sqrt( Core::K_BOLTZMANN * temperature_K / (collisionGasMass_kg_) );
        Core::Vector velocityBgMolecule = { (mole.getComPos().x()-bgMole.getComPos().x())/fabs(mole.getComPos().x()-bgMole.getComPos().x()) * 
                                            rndSource->normalRealRndValue() * vrStdevBgMolecule - particle.getVelocity().x(),
                                            0,
                                            0};
        bgMole.setComVel(velocityBgMolecule);

        // // calculate random point on sphere
        // // as follows:
        // // draw random number in as long until magnitude is less than 1
        // // normalize result
        // double circleVectorMagnitude = 0;
        // Core::Vector circleVector(0,0,0);
        // double directionAngle = 0;
        // do{
        //     circleVector = Core::Vector{
        //           (rndSource->uniformRealRndValue() * 2 - 1)
        //         , (rndSource->uniformRealRndValue() * 2 - 1)
        //         , (rndSource->uniformRealRndValue() * 2 - 1)
        //     };
        //     circleVectorMagnitude = circleVector.magnitude();
        //     circleVector = spawnRad / circleVectorMagnitude * circleVector;
        //     directionAngle  = std::acos(
        //         ( (-1. * circleVector) * velocityBgMolecule)
        //         / ( spawnRad * velocityBgMolecule.magnitude() )
        //     );
        // }while(circleVectorMagnitude > 1 || directionAngle >  angleThetaScaling_ * collisionTheta);


        // NO ROTATION

        // // rotate it randomly
        // bgMole.setAngles(Core::Vector(rndSource->uniformRealRndValue(),
        //                             rndSource->uniformRealRndValue(),
        //                             rndSource->uniformRealRndValue()));

        // // Give molecule a random orientation:
        // mole.setAngles(Core::Vector(rndSource->uniformRealRndValue(),
        //                             rndSource->uniformRealRndValue(),
        //                             rndSource->uniformRealRndValue()));

        std::vector<CollisionModel::Molecule*> moleculesPtr = {&mole, &bgMole};

        // possible check for energy conservation
        std::vector<Core::Vector> startVelocity;
        double startEnergy = 0;
        for(auto* molecule : moleculesPtr){
            startVelocity.push_back(molecule->getComVel());
            startEnergy += 0.5 * molecule->getMass() * molecule->getComVel().magnitudeSquared();
        }

        // Call the sub-integrator
        double finalTime = integrationTime_; //  final integration time in seconds
        double timeStep = subTimeStep_; // step size in seconds

        // Core::Vector bgMole_originalComPos = bgMole.getComPos();
        // Core::Vector bgMole_originalComVel = bgMole.getComVel();
        // Core::Vector mole_originalComPos = mole.getComPos();
        // Core::Vector mole_originalComVel = mole.getComVel();

        // trajectorySuccess = leapfrogIntern(moleculesPtr, timeStep, finalTime, collisionRadius);
        trajectorySuccess = rk4InternAdaptiveStep(moleculesPtr, timeStep, finalTime, collisionRadius);
       
        double endEnergy = 0;
        for(auto* molecule : moleculesPtr){
            endEnergy += 0.5 * molecule->getMass() * molecule->getComVel().magnitudeSquared();
        }
        std::cout << startEnergy << " " << endEnergy << std::endl;
        if(endEnergy*0.90 >= startEnergy){
            trajectorySuccess = false;
            dt = dt*0.98;
        }
        if(trajectorySuccess){
            particle.setVelocity(mole.getComVel() + particle.getVelocity());
        }
        ++iterations;
    }while(!trajectorySuccess && iterations < 1);

    if(trajectorySuccess == false){
        std::cerr << "No trajectory that hit the collision sphere was found.\n";
    }
}

void CollisionModel::MDInteractionsModelPreconstructed::modifyPosition(Core::Vector& /*position*/, Core::Particle& /*particle*/, double /*dt*/) {

}

bool CollisionModel::MDInteractionsModelPreconstructed::leapfrogIntern(std::vector<CollisionModel::Molecule*> moleculesPtr, double dt, double finalTime, double requiredRad){

    // static std::ofstream myfile;
    // if(!myfile.is_open()){
    //     myfile.open("exampleLF.txt");
    // }
    bool wasHit = false;

    // distances need to be saved so integration can be stopped if particles leave 
    // the domain of interest 
    std::vector<double> startDistances;
    size_t moleculesPtr_size = moleculesPtr.size();
    for(size_t i = 0; i < moleculesPtr_size; ++i){
        for(size_t j = i+1; j < moleculesPtr_size; ++j){
            startDistances.push_back((moleculesPtr.at(i)->getComPos() - moleculesPtr.at(j)->getComPos()).magnitude());
        }
    }


    int nSteps = int(round(finalTime/dt));

    std::vector<Core::Vector> forceMolecules(moleculesPtr_size);
    forceFieldMD(moleculesPtr, forceMolecules);

    // do the first half step for the velocity, as per leapfrog definition
    double energyStart = 0;
    size_t i = 0;
    for(auto* molecule : moleculesPtr){
        energyStart += 0.5 * molecule->getComVel().magnitudeSquared() * molecule->getMass();
        Core::Vector newComVel =  molecule->getComVel() + forceMolecules.at(i) / molecule->getMass() * dt/2;
        molecule->setComVel(newComVel);
        i++;
    }

    // start the actual leapfrog iteration
    for (int j = 0; j < nSteps; j++){
        // time step for the new position
        i = 0;
        double energyEnd = 0;
        for(auto* molecule : moleculesPtr){
            Core::Vector newComPos =  molecule->getComPos() + molecule->getComVel() * dt;
            molecule->setComPos(newComPos);
            energyEnd += 0.5 * molecule->getComVel().magnitudeSquared() * molecule->getMass();
            i++;
        }
        size_t index = 0;
        for(size_t k = 0; k < moleculesPtr_size; ++k){
            for(size_t l = k+1; l < moleculesPtr_size; ++l){
                if((moleculesPtr.at(l)->getComPos() - moleculesPtr.at(k)->getComPos()).magnitude() > startDistances.at(index++)){
                    return wasHit;
                }
                else if((moleculesPtr.at(l)->getComPos() - moleculesPtr.at(k)->getComPos()).magnitude() <= requiredRad){
                    wasHit=true;
                }
            }
        }

        // recalculate the force
        forceFieldMD(moleculesPtr, forceMolecules);
        i = 0;
        // time step for the new velocity
        for(auto* molecule : moleculesPtr){
            Core::Vector newComVel =  molecule->getComVel() + forceMolecules.at(i) / molecule->getMass() * dt;
            molecule->setComVel(newComVel);
            i++;
        }
        // i = 0;
        // for(auto* molecule : moleculesPtr){
        //     myfile << molecule->getComVel() << " ";
        //     i++;
        // }
        // myfile << std::endl;
    }
    return false;

}

bool CollisionModel::MDInteractionsModelPreconstructed::rk4Intern(std::vector<CollisionModel::Molecule*> moleculesPtr, double dt, double finalTime,
                                                                    double requiredRad){


    int nSteps = int(round(finalTime/dt));
    size_t nMolecules = moleculesPtr.size();
    std::vector<Core::Vector> forceMolecules(nMolecules);

    bool wasHit = false;
    std::vector<double> startDistances;
    for(size_t i = 0; i < nMolecules; ++i){
        for(size_t j = i+1; j < nMolecules; ++j){
            startDistances.push_back((moleculesPtr[i]->getComPos() - moleculesPtr[j]->getComPos()).magnitude());
        }
    }

    size_t i = 0;

    for (int j = 0; j < nSteps; j++){

        std::vector<Core::Vector> velocityMolecules(nMolecules);
        std::vector<Core::Vector> positionMolecules(nMolecules);
        i = 0;
        for(auto* molecule : moleculesPtr){
            velocityMolecules.at(i) = molecule->getComVel();
            positionMolecules.at(i) = molecule->getComPos();
            i++;
        }
        
        std::vector<Core::Vector> initialPositionMolecules(nMolecules);
        std::vector<Core::Vector> initialVelocityMolecules(nMolecules);
        for(size_t k = 0; k < nMolecules; k++){
            initialPositionMolecules.at(k) = Core::Vector( positionMolecules.at(k).x(), positionMolecules.at(k).y(),positionMolecules.at(k).z() );
            initialVelocityMolecules.at(k) = Core::Vector( velocityMolecules.at(k).x(), velocityMolecules.at(k).y(), velocityMolecules.at(k).z() );
        }
        
        double length[3] = {1./2, 1./2, 1};
        double mass[nMolecules];
        i = 0;
        for(auto* molecule : moleculesPtr){
            mass[i] = molecule->getMass();
            i++;
        }
        forceFieldMD(moleculesPtr, forceMolecules);

        std::array<std::array<Core::Vector, 2>, 4> k;
        std::array<std::array<Core::Vector, 2>, 4> l;


        for(size_t q = 0; q < nMolecules; q++){
            k[0][q] = forceMolecules.at(q) * dt / mass[q];
            l[0][q] = velocityMolecules.at(q) * dt;
        }

        for(size_t n = 1; n < 4; n++){
            i = 0;
            for(auto* molecule : moleculesPtr){
                positionMolecules.at(i) = initialPositionMolecules.at(i) + l[n-1][i]*length[i-1];
                molecule->setComPos(positionMolecules.at(i));
                i++;
            }

            forceFieldMD(moleculesPtr, forceMolecules);

            for(i = 0; i < nMolecules; i++){
                k[n][i] = forceMolecules.at(i) * dt / mass[i];
                l[n][i] = (velocityMolecules.at(i) + k[n-1][i]*length[n-1])*dt;
            }

        }

        i = 0;
        for(auto* molecule : moleculesPtr){
            Core::Vector newComPos = initialPositionMolecules.at(i) + (l[0][i]+ l[1][i]*2 + l[2][i]*2 + l[3][i]) * 1./6;
            molecule->setComPos(newComPos);
            Core::Vector newComVel = initialVelocityMolecules.at(i) + (k[0][i]+ k[1][i]*2 + k[2][i]*2 + k[3][i]) * 1./6;
            molecule->setComVel(newComVel);
            i++;
        }

        size_t index = 0;
        for(size_t k = 0; k < nMolecules; ++k){
            for(size_t l = k+1; l < nMolecules; ++l){
                if((moleculesPtr[l]->getComPos() - moleculesPtr[k]->getComPos()).magnitude() > startDistances[index++]){
                    return wasHit;
                }
                if((moleculesPtr[l]->getComPos() - moleculesPtr[k]->getComPos()).magnitude() <= requiredRad){
                    wasHit=true;
                }
            }
        }
    }

    return false;
}


bool CollisionModel::MDInteractionsModelPreconstructed::rk4InternAdaptiveStep(std::vector<CollisionModel::Molecule*> moleculesPtr, double dt, double finalTime,
                                                                    double requiredRad){

    double integrationTimeSum = 0;
    size_t nMolecules = moleculesPtr.size();
    std::vector<Core::Vector> forceMolecules(nMolecules);

    size_t i = 0;
    int steps = 0;
    double distance = 0.0;

    // distances need to be saved so integration can be stopped if particles leave 
    // the domain of interest 
    bool wasHit = false;
    std::vector<double> startDistances;
    for(size_t i = 0; i < nMolecules; ++i){
        for(size_t j = i+1; j < nMolecules; ++j){
            startDistances.push_back((moleculesPtr[i]->getComPos() - moleculesPtr[j]->getComPos()).magnitude());
        }
    }
    std::vector<Core::Vector> velocityMolecules(nMolecules);
    std::vector<Core::Vector> positionMolecules(nMolecules);
    std::vector<Core::Vector> initialPositionMolecules(nMolecules);
    std::vector<Core::Vector> initialVelocityMolecules(nMolecules);

    double weight[5][6] = { 
                            {1./4, 0, 0, 0, 0, 0},
                            {3./32, 9./32, 0, 0, 0, 0},
                            {1932./2197, -7200./2197, 7296./2197, 0, 0, 0},
                            {439./216, -8, 3680./513, -845./4104, 0, 0},
                            {-8./27, 2, -3544./2565, 1859./4104, -11./40, 0}};
    double mass[2];
    std::array<std::array<Core::Vector, 2>, 6> k;
    std::array<std::array<Core::Vector, 2>, 6> l;

    while(integrationTimeSum < finalTime){
        
        i = 0;
        for(auto* molecule : moleculesPtr){
            velocityMolecules[i] = molecule->getComVel();
            positionMolecules[i] = molecule->getComPos();
            i++;
        }


        for(size_t k = 0; k < nMolecules; k++){
            initialPositionMolecules[k] = Core::Vector( positionMolecules[k].x(), positionMolecules[k].y(),positionMolecules[k].z() );
            initialVelocityMolecules[k] = Core::Vector( velocityMolecules[k].x(), velocityMolecules[k].y(), velocityMolecules[k].z() );
        }

        i = 0;
        for(auto* molecule : moleculesPtr){
            mass[i] = molecule->getMass();
            i++;
        }

        forceFieldMD(moleculesPtr, forceMolecules);

        for(size_t q = 0; q < nMolecules; q++){
            k[0][q] = forceMolecules[q] * dt / mass[q];
            l[0][q] = velocityMolecules[q] * dt;
        }

        for(size_t n = 1; n < 6; n++){
            i = 0;
            for(auto* molecule : moleculesPtr){
                positionMolecules[i] = initialPositionMolecules[i];
                // molecule->setComPos(positionMolecules[i]);
                i++;
            }
            for(size_t m = 0; m < 6; m++){
                i = 0;
                for(auto* molecule : moleculesPtr){
                    positionMolecules[i] += l[m][i]*weight[n-1][m];
                    molecule->setComPos(positionMolecules[i]);
                    i++;
                }
            }

            forceFieldMD(moleculesPtr, forceMolecules);
            
            
            for(i = 0; i < nMolecules; i++){
                k[n][i] = forceMolecules[i] * dt / mass[i];
                l[n][i] = velocityMolecules[i];
                for(size_t m = 0; m < 6; m++){
                    l[n][i] += k[m][i]*weight[n-1][m];
                }
                l[n][i] = l[n][i]*dt;
            }

        }

        for(size_t k = 0; k < nMolecules; ++k){
            for(size_t l = k+1; l < nMolecules; ++l){
                distance = (moleculesPtr[l]->getComPos() - moleculesPtr[k]->getComPos()).magnitude();
            }
        }

        i = 0;
        std::array<Core::Vector, 2> newComVelOrder5;
        std::array<Core::Vector, 2> newComPosOrder4; 
        std::array<Core::Vector, 2> newComVelOrder4;
        for(i = 0; i < 2; i++){
            newComVelOrder5[i] = initialVelocityMolecules[i] + (k[0][i] * 16./135 + k[2][i] * 6656./12825 + k[3][i] * 28561./56430 + k[4][i] * (-9./50) + k[5][i] * 2./55);
            newComPosOrder4[i] = initialPositionMolecules[i] + (l[0][i] * 25./216 + l[2][i] * 1408./2565 + l[3][i] * 2197./4104 + l[4][i] * (-1./5));
            newComVelOrder4[i] = initialVelocityMolecules[i] + (k[0][i] * 25./216 + k[2][i] * 1408./2565 + k[3][i] * 2197./4104 + k[4][i] * (-1./5));
        }

        std::array<double,2> RX;
        std::array<double,2> RY;
        std::array<double,2> RZ;
        for(size_t p = 0; p < 2; p++){
            if(fabs(newComVelOrder5[p].x()) != 0)
                RX[p] = fabs(newComVelOrder4[p].x() - newComVelOrder5[p].x()) / fabs(newComVelOrder5[p].x());
            else
                RX[p] = 0;
            if(fabs(newComVelOrder5[p].y()) != 0)
                RY[p] = fabs(newComVelOrder4[p].y() - newComVelOrder5[p].y()) / fabs(newComVelOrder5[p].y());
            else
                RY[p] = 0;
            if(fabs(newComVelOrder5[p].z()) != 0)
                RZ[p] = fabs(newComVelOrder4[p].z() - newComVelOrder5[p].z()) / fabs(newComVelOrder5[p].z());
            else
                RZ[p] = 0;            
        }

        double globalX = std::max({RX[0],RX[1]});
        double globalY = std::max({RY[0],RY[1]});
        double globalZ = std::max({RZ[0],RZ[1]});
        
        double tolerance = 1e-8;
        double globalR = std::max({globalX, globalY, globalZ});
        if (globalR == 0){
            globalR = 1e-15;
        }
        
        double globalDelta = 0.84 * std::pow((tolerance/globalR), 1./4);
        integrationTimeSum += dt;
        i = 0;
        for(auto* molecule : moleculesPtr){
            if(trajectoryRecordingActive_ == true && molecule->getMolecularStructureName() == collisionMolecule_){
                writeTrajectory(distance, molecule->getComPos(), molecule->getComVel(),forceMolecules, false, trajectoryOutputStream_.get(), integrationTimeSum, dt);
            }
            molecule->setComPos(newComPosOrder4[i]);
            molecule->setComVel(newComVelOrder4[i]);
            i++;
        }
        steps++;
        dt = dt * globalDelta;

        size_t index = 0;
        for(size_t k = 0; k < nMolecules; ++k){
            for(size_t l = k+1; l < nMolecules; ++l){
                if((moleculesPtr[l]->getComPos() - moleculesPtr[k]->getComPos()).magnitude() > startDistances[index++]){
                    if(trajectoryRecordingActive_ == true && moleculesPtr[l]->getMolecularStructureName() == collisionMolecule_ && wasHit == true){
                        writeTrajectory((moleculesPtr[l]->getComPos() - moleculesPtr[k]->getComPos()).magnitude(),
                                        moleculesPtr[l]->getComPos(), moleculesPtr[l]->getComVel(), forceMolecules, true, trajectoryOutputStream_.get(), integrationTimeSum, dt);
                    }
                    return wasHit;
                }
                if((moleculesPtr[l]->getComPos() - moleculesPtr[k]->getComPos()).magnitude() <= requiredRad){
                    wasHit=true;
                }
            }
        }
    }
  
    return false;
}


void CollisionModel::MDInteractionsModelPreconstructed::forceFieldMD(std::vector<CollisionModel::Molecule*>& moleculesPtr, std::vector<Core::Vector>& forceMolecules){

    //std::vector<Core::Vector> forceMolecules(nMolecules); 
    // save all the forces acting on each molecule
    CollisionModel::Molecule* ion = moleculesPtr[0];
    CollisionModel::Molecule* bgGas = moleculesPtr[1];
    forceMolecules[0] = Core::Vector(0.0, 0.0, 0.0);
    forceMolecules[1] = Core::Vector(0.0, 0.0, 0.0);
    bool isN2 = false;

    if(bgGas->getMolecularStructureName()=="N2"){
        isN2 = true;
    }

    // construct E-field acting on the molecule
    std::array<double, 3> eField = {0., 0., 0.};
    std::array<double, 6> eFieldDerivative = {0., 0., 0., 0., 0., 0.};


    /* 
    * therefore we need the interaction between each atom of a molecule with the atoms of the
    * other one 
    */ 
    for(auto& atomI : ion->getAtoms()){
        for(auto& atomJ : bgGas->getAtoms()){
            
            // First contribution: Lennard-Jones potential 
            // This always contributes to the experienced force 
            Core::Vector absPosAtomI = ion->getComPos() + atomI->getRelativePosition();
            Core::Vector absPosAtomJ = bgGas->getComPos() + atomJ->getRelativePosition();


            Core::Vector distance = absPosAtomI - absPosAtomJ;
            // avoid singularity at zero distance 
            if(distance.magnitude() < 1E-25){
                forceMolecules[0] += Core::Vector(1e-10, 1e-10, 1e-10);
                forceMolecules[1] += Core::Vector(1e-10, 1e-10, 1e-10) * (-1);
                break;
            }
            // cut-off
            if(distance.magnitude() > 1E20){
                return;
            }
            double distanceSquared = distance.magnitudeSquared();
            double distanceSquaredInverse = 1./distanceSquared;
            double sigma = CollisionModel::Atom::calcLJSig(*atomI, *atomJ);
            double sigma6 = sigma * sigma * sigma * sigma * sigma * sigma;
            double epsilon = CollisionModel::Atom::calcLJEps(*atomI, *atomJ);
            double ljFactor = 24 * epsilon * distanceSquaredInverse*distanceSquaredInverse*distanceSquaredInverse*distanceSquaredInverse * 
                                (2 * distanceSquaredInverse*distanceSquaredInverse*distanceSquaredInverse * sigma6 * sigma6 - sigma6);
            // calculate the force that acts on the atoms and add it to the overall force on the molecule
            Core::Vector atomForce;
            atomForce.x(distance.x() * ljFactor);
            atomForce.y(distance.y() * ljFactor);
            atomForce.z(distance.z() * ljFactor);
            forceMolecules[0] += atomForce;
            forceMolecules[1] += atomForce * (-1);

            // Second contribution: C4 ion-induced dipole potential
            // This requires an ion and one neutrally charged molecule to be present
            double distanceCubed = distanceSquared * sqrt(distanceSquared);
            double currentCharge = 0;
            // Check if one of the molecules is an ion and the other one is not
            if(int(atomI->getCharge()/Core::ELEMENTARY_CHARGE) != 0 && 
                moleculesPtr[1]->getIsIon() == false &&
                moleculesPtr[1]->getIsDipole() == false){
                currentCharge = atomI->getCharge();

            }else if (moleculesPtr[0]->getIsIon() == false && 
                        int(atomJ->getCharge()/Core::ELEMENTARY_CHARGE) != 0 &&
                        moleculesPtr[0]->getIsDipole() == false){
                currentCharge = atomJ->getCharge();
            }
            
            if(distance.magnitude() <= 30e-10){
                eField[0] += distance.x() * currentCharge / distanceCubed; // E-field in x
                eField[1] += distance.y() * currentCharge / distanceCubed; // E-field in y
                eField[2] += distance.z() * currentCharge / distanceCubed; // E-field in z
                
                // derivative x to x
                eFieldDerivative[0] += currentCharge / distanceCubed - 
                                        3 * currentCharge * distance.x() * distance.x() / (distanceCubed * distanceSquared); 
                // derivative x to y
                eFieldDerivative[1] += -3 * currentCharge * distance.x() * distance.y() / (distanceCubed * distanceSquared);
                // derivative y to y
                eFieldDerivative[2] += currentCharge / distanceCubed - 
                                        3 * currentCharge * distance.y() * distance.y() / (distanceCubed * distanceSquared);
                // derivative y to z
                eFieldDerivative[3] += -3 * currentCharge * distance.y() * distance.z() / (distanceCubed * distanceSquared);
                // derivative z to z
                eFieldDerivative[4] += currentCharge / distanceCubed - 
                                        3 * currentCharge * distance.z() * distance.z() / (distanceCubed * distanceSquared);
                // derivative x to z
                eFieldDerivative[5] += -3 * currentCharge * distance.x() * distance.z() / (distanceCubed * distanceSquared);
            }
            

            // Third contribution: ion <-> permanent dipole potential
            // This requires an ion and a dipole to be present 
            double dipoleDistanceScalar = 0;
            double dipoleX = 0, dipoleY = 0, dipoleZ = 0;
            currentCharge = 0;
            if(int(atomI->getCharge()/Core::ELEMENTARY_CHARGE) != 0 && 
                moleculesPtr[1]->getIsDipole() == true){

                currentCharge = atomI->getCharge();
                dipoleX = moleculesPtr[1]->getDipole().x();
                dipoleY = moleculesPtr[1]->getDipole().y();
                dipoleZ = moleculesPtr[1]->getDipole().z();
                dipoleDistanceScalar =  dipoleX * distance.x() + 
                                        dipoleY * distance.y() + 
                                        dipoleZ * distance.z();

            }else if (moleculesPtr[0]->getIsDipole() == true && 
                        int(atomJ->getCharge()/Core::ELEMENTARY_CHARGE) != 0){

                currentCharge = atomJ->getCharge();
                dipoleX = moleculesPtr[0]->getDipole().x();
                dipoleY = moleculesPtr[0]->getDipole().y();
                dipoleZ = moleculesPtr[0]->getDipole().z();
                dipoleDistanceScalar =  dipoleX * distance.x() + 
                                        dipoleY * distance.y() + 
                                        dipoleZ * distance.z();
            }
            // Core::Vector ionDipoleForce;
            // ionDipoleForce.x(-currentCharge * 1./Core::ELECTRIC_CONSTANT * 
            //                     (1./distanceCubed * dipoleX - 
            //                     3 * dipoleDistanceScalar * 1./(distanceCubed*distanceSquared) * distance.x()) );
            // ionDipoleForce.y(-currentCharge * 1./Core::ELECTRIC_CONSTANT * 
            //                     (1./distanceCubed * dipoleY - 
            //                     3 * dipoleDistanceScalar * 1./(distanceCubed*distanceSquared) * distance.y()) );
            // ionDipoleForce.z(-currentCharge * 1./Core::ELECTRIC_CONSTANT * 
            //                     (1./distanceCubed * dipoleZ - 
            //                     3 * dipoleDistanceScalar * 1./(distanceCubed*distanceSquared) * distance.z()) );
            // forceMolecules[0] += ionDipoleForce;
            // forceMolecules[1] += ionDipoleForce * (-1);

            // Fourth contribution: quadrupole moment if background gas is N2 
            // This requires an ion and N2 to be present 
            double partialChargeN2 = 0;
            if(int(atomI->getCharge()/Core::ELEMENTARY_CHARGE) != 0 && 
                isN2 == true){

                currentCharge = atomI->getCharge();
                partialChargeN2 = atomJ->getPartCharge();

            }else if (isN2 == true && 
                        int(atomJ->getCharge()/Core::ELEMENTARY_CHARGE) != 0){

                currentCharge = atomJ->getCharge();
                partialChargeN2 = atomI->getPartCharge();
        
            }
            // Core::Vector quadrupoleForce;
            // quadrupoleForce.x(-currentCharge * partialChargeN2 * 1./Core::ELECTRIC_CONSTANT * distance.x() / distanceCubed );
            // quadrupoleForce.y(-currentCharge * partialChargeN2 * 1./Core::ELECTRIC_CONSTANT * distance.y() / distanceCubed );
            // quadrupoleForce.z(-currentCharge * partialChargeN2 * 1./Core::ELECTRIC_CONSTANT * distance.z() / distanceCubed );
            // forceMolecules[0] += quadrupoleForce;
            // forceMolecules[1] += quadrupoleForce * (-1);

        }
    }

    // add the C4 ion-induced dipole force contribution 
    Core::Vector ionInducedForce;
    ionInducedForce.x(1./Core::ELECTRIC_CONSTANT * collisionGasPolarizability_m3_ * 
                        (eField[0]*eFieldDerivative[0] + eField[1]*eFieldDerivative[1] + eField[2]*eFieldDerivative[5]));
    ionInducedForce.y(1./Core::ELECTRIC_CONSTANT * collisionGasPolarizability_m3_ * 
                        (eField[0]*eFieldDerivative[1] + eField[1]*eFieldDerivative[2] + eField[2]*eFieldDerivative[3]));
    ionInducedForce.z(1./Core::ELECTRIC_CONSTANT * collisionGasPolarizability_m3_ * 
                        (eField[0]*eFieldDerivative[5] + eField[1]*eFieldDerivative[3] + eField[2]*eFieldDerivative[4]));
    forceMolecules[0] += ionInducedForce;
    forceMolecules[1] += ionInducedForce * (-1);
    
}