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
                                                        std::unique_ptr<CollisionModel::AbstractMDForceField> forceField,
                                                        std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection) :
        MDInteractionsModel(
        getConstantDoubleFunction(staticPressure),
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
                                                        std::unique_ptr<CollisionModel::AbstractMDForceField> forceField,
                                                        std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection) :

        MDInteractionsModel(
                std::move(pressureFunction),
                std::move(velocityFunction),
                getConstantDoubleFunction(staticTemperature),
                collisionGasMassAmu,
                collisionGasDiameterM,
                collisionMolecule,
                integrationTime,
                subTimeStep,
                collisionRadiusScaling,
                angleThetaScaling,
                spawnRadius,
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
        pressureFunction_(std::move(pressureFunction)),
        velocityFunction_(std::move(velocityFunction)),
        temperatureFunction_(std::move(temperatureFunction)),
        collisionGasMass_kg_(collisionGasMassAmu*Core::AMU_TO_KG),
        collisionGasDiameter_m_(collisionGasDiameterM),
        collisionMolecule_(collisionMolecule),
        integrationTime_(integrationTime),
        subTimeStep_(subTimeStep),
        collisionRadiusScaling_(collisionRadiusScaling),
        angleThetaScaling_(angleThetaScaling),
        spawnRadius_(spawnRadius),
        forceField_(std::move(forceField)),
        molecularStructureCollection_(std::move(molecularStructureCollection)) { }

/**
 * Activates trajectory writing and sets trajectory writer configuration
 * @param trajectoryFileName Trajectory Output filename
 * @param trajectoryDistance Distance between ion and background gas in m after which trajectory gets recorded
 */
void CollisionModel::MDInteractionsModel::setTrajectoryWriter(const std::string& trajectoryFileName,
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
/**
 * Returns sign of a number (+1, -1) or 0.
*/
double CollisionModel::MDInteractionsModel::calcSign(double value){
    if(value > 0){
        return 1.;
    }else if(value < 0){
        return -1.;
    }else{
        return 0;
    }
}

/**
 * Writes trajectory data to a predefined file. Individual collisions are separated by a line containing '###'. 
 * 
 * Ouput includes: position of background gas, distance between the two molecules, 
 * integration time, velocity of the background gas, force acting on the background gas and 
 * timestep length
 */
void CollisionModel::MDInteractionsModel::writeTrajectory(double distance, Core::Vector positionBgMolecule, Core::Vector velocityBgMolecule, 
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

void CollisionModel::MDInteractionsModel::initializeModelParticleParameters(Core::Particle& /*ion*/) const {

}

void CollisionModel::MDInteractionsModel::updateModelParticleParameters(Core::Particle& /*ion*/) const {

}

/**
 * Updates trajectory recording if timestep recording parameter is exceeded
*/
void CollisionModel::MDInteractionsModel::updateModelTimestepParameters(int timestep, double /*time*/) {
    if (modelRecordsTrajectories_ && timestep >= recordTrajectoryStartTimeStep_){
        trajectoryRecordingActive_ = true;
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

        double pi = 3.1415; 
        // // rotate it randomly
        // bgMole.setAngles(Core::Vector(rndSource->uniformRealRndValue()*2*pi-pi,
        //                               rndSource->uniformRealRndValue()*2*pi-pi,
        //                               rndSource->uniformRealRndValue()*2*pi-pi));

        // Give molecule a random orientation:
        mole.setAngles(Core::Vector(rndSource->uniformRealRndValue()*2*pi-pi,
                                    rndSource->uniformRealRndValue()*2*pi-pi,
                                    rndSource->uniformRealRndValue()*2*pi-pi));

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

        trajectorySuccess = rk4InternAdaptiveStep(moleculesPtr, timeStep, finalTime, collisionRadius, tolerance);
        // trajectorySuccess = leapfrogIntern(moleculesPtr, timeStep, finalTime, collisionRadius);
       
        double endEnergy = 0;
        for(auto* molecule : moleculesPtr){
            endEnergy += 0.5 * molecule->getMass() * molecule->getComVel().magnitudeSquared();
        }
        // check if energy is conserved up to 10% 
        // if not halve the starting timestep length 
        if(endEnergy*0.90 >= startEnergy){
            std::cout << "Energy not conserved: " << startEnergy << " " << endEnergy << std::endl;
            trajectorySuccess = false;
            dt = dt/2;
            // tolerance /= 2;
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

void CollisionModel::MDInteractionsModel::modifyPosition(Core::Vector& /*position*/, Core::Particle& /*particle*/, double /*dt*/) {

}

/**
 * Leapfrog method to integrate trajectories of particles involved in a collision. 
 * The leapfrog method is of second order and symplectic. 
 * @param moleculesPtr collection of molecule pointer  
 * @param dt timestep length 
 * @param finalTime maximum integration time 
 * @param requiredRad radius defining the collision sphere, i.e. the distance that needs to be undercut for 
 * a collision to be considered (same radius which is used to estimate the collision probability)
 */
bool CollisionModel::MDInteractionsModel::leapfrogIntern(std::vector<CollisionModel::Molecule*> moleculesPtr, double dt, double finalTime, double requiredRad){

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
    forceField_->calculateForceField(moleculesPtr, forceMolecules);

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
        forceField_->calculateForceField(moleculesPtr, forceMolecules);
        i = 0;
        // time step for the new velocity
        for(auto* molecule : moleculesPtr){
            Core::Vector newComVel =  molecule->getComVel() + forceMolecules.at(i) / molecule->getMass() * dt;
            molecule->setComVel(newComVel);
            i++;
        }
    }
    return false;

}

/**
 * Runge-Kutta 4 method to integrate trajectories of particles involved in a collision.
 * The RK4 is of fourth order.
 * This integrator should NOT be used except for testing purposes as the adaptive step size method is faster.
 * @param moleculesPtr collection of molecule pointer 
 * @param dt timestep length 
 * @param finalTime maximum integration time 
 * @param requiredRad radius defining the collision sphere, i.e. the distance that needs to be undercut for 
 * a collision to be considered (same radius which is used to estimate the collision probability)
 */
bool CollisionModel::MDInteractionsModel::rk4Intern(std::vector<CollisionModel::Molecule*> moleculesPtr, double dt, double finalTime,
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
        double mass[2];
        i = 0;
        for(auto* molecule : moleculesPtr){
            mass[i] = molecule->getMass();
            i++;
        }
        forceField_->calculateForceField(moleculesPtr, forceMolecules);

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

            forceField_->calculateForceField(moleculesPtr, forceMolecules);

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

/**
 * Adaptive step size Runge-Kutta-Fehlberg 45 method to integrate trajectories of particles involved in a collision. 
 * This method is of fourth order and uses error control of fifth order on the velocity to adaptively 
 * increase or decrease the timestep length reducing the overall computation time. 
 * @param moleculesPtr collection of molecule pointer 
 * @param dt timestep length 
 * @param finalTime maximum integration time 
 * @param requiredRad radius defining the collision sphere, i.e. the distance that needs to be undercut for 
 * a collision to be considered (same radius which is used to estimate the collision probability)
 * @param tolerance defines the allowed error threshold  to control the timestep lengths 
 */
bool CollisionModel::MDInteractionsModel::rk4InternAdaptiveStep(std::vector<CollisionModel::Molecule*> moleculesPtr, double dt, double finalTime,
                                                                    double requiredRad, double tolerance){

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
    for(i = 0; i < nMolecules; ++i){
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
    std::array<Core::Vector, 2> newComVelOrder5;
    std::array<Core::Vector, 2> newComPosOrder4; 
    std::array<Core::Vector, 2> newComVelOrder4;
    std::array<double, 2> R;
    double globalR, globalDelta;
    // double tolerance = 1e-8;
    double pi = 3.14159;
    Core::RandomSource* rndSource = Core::globalRandomGeneratorPool->getThreadRandomSource();
    Core::Vector nitrogenOne;
    Core::Vector nitrogenTwo; 
    Core::Vector nitrogenAngles = {rndSource->uniformRealRndValue()*2*pi-pi, 
                                rndSource->uniformRealRndValue()*2*pi-pi, 
                                rndSource->uniformRealRndValue()*2*pi-pi};
    double I = 0;
    double angularVelocity = 0;
    if(moleculesPtr[1]->getMolecularStructureName()=="N2"){
        nitrogenOne = molecularStructureCollection_.at(moleculesPtr[1]->getMolecularStructureName())->getAtoms().at(0)->getRelativePosition();
        nitrogenTwo = molecularStructureCollection_.at(moleculesPtr[1]->getMolecularStructureName())->getAtoms().at(1)->getRelativePosition();
        I = CollisionModel::MolecularStructure::getMomentOfInertia(nitrogenOne.x(), nitrogenTwo.x(), 
                                                                    moleculesPtr[1]->getMass()/2, moleculesPtr[1]->getMass()/2);
        angularVelocity = CollisionModel::MolecularStructure::getAngularVelocity(temperatureFunction_(moleculesPtr[1]->getComPos()), I);
    }
    moleculesPtr[1]->setAngles(nitrogenAngles);

    while(integrationTimeSum < finalTime){
        
        i = 0;
        for(auto* molecule : moleculesPtr){
            velocityMolecules[i] = molecule->getComVel();
            positionMolecules[i] = molecule->getComPos();
            initialPositionMolecules[i] = molecule->getComPos();
            initialVelocityMolecules[i] = molecule->getComVel();
            mass[i] = molecule->getMass();
            i++;
        }

        forceField_->calculateForceField(moleculesPtr, forceMolecules);

        for(size_t q = 0; q < nMolecules; q++){
            k[0][q] = forceMolecules[q] * dt / mass[q];
            l[0][q] = velocityMolecules[q] * dt;
        }

        for(size_t n = 1; n < 6; n++){
            for(i = 0; i < nMolecules; i++){
                positionMolecules[i] = initialPositionMolecules[i];
                
            }
            for(size_t m = 0; m < 6; m++){
                i = 0;
                for(auto* molecule : moleculesPtr){
                    positionMolecules[i] += l[m][i]*weight[n-1][m];
                    molecule->setComPos(positionMolecules[i]);
                    i++;
                }
            }

            forceField_->calculateForceField(moleculesPtr, forceMolecules);
            
            for(i = 0; i < nMolecules; i++){
                k[n][i] = forceMolecules[i] * dt / mass[i];
                l[n][i] = velocityMolecules[i];
                for(size_t m = 0; m < 6; m++){
                    l[n][i] += k[m][i]*weight[n-1][m];
                }
                l[n][i] = l[n][i]*dt;
            }

        }

        for(size_t b = 0; b < nMolecules; ++b){
            for(size_t z = b+1; z < nMolecules; ++z){
                distance = (moleculesPtr[z]->getComPos() - moleculesPtr[b]->getComPos()).magnitude();
            }
        }

        
        for(i = 0; i < 2; i++){
            newComVelOrder5[i] = initialVelocityMolecules[i] + (k[0][i] * 16./135 + k[2][i] * 6656./12825 + k[3][i] * 28561./56430 + k[4][i] * (-9./50) + k[5][i] * 2./55);
            newComPosOrder4[i] = initialPositionMolecules[i] + (l[0][i] * 25./216 + l[2][i] * 1408./2565 + l[3][i] * 2197./4104 + l[4][i] * (-1./5));
            newComVelOrder4[i] = initialVelocityMolecules[i] + (k[0][i] * 25./216 + k[2][i] * 1408./2565 + k[3][i] * 2197./4104 + k[4][i] * (-1./5));
        }

        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wfloat-equal"
        for(size_t p = 0; p < 2; p++){
            if(fabs(newComVelOrder5[p].magnitude()) != 0)
                R[p] = fabs(newComVelOrder4[p].magnitude()-newComVelOrder5[p].magnitude())/fabs(newComVelOrder5[p].magnitude());
            else
                R[p] = 0;
                   
        }

        globalR = std::max({R[0],R[1]});

        if (globalR == 0){
            globalR = 1e-15;
        }
        #pragma GCC diagnostic pop

        globalDelta = 0.84 * std::pow((tolerance/globalR), 1./4);
        integrationTimeSum += dt;
        i = 0;
        for(auto* molecule : moleculesPtr){
            if(trajectoryRecordingActive_ == true && molecule->getMolecularStructureName() == collisionMolecule_){
                writeTrajectory(distance, molecule->getComPos(), molecule->getComVel(),forceMolecules, false, trajectoryOutputStream_.get(), integrationTimeSum, dt);
            }
            
            molecule->setComPos(newComPosOrder4[i]);
            molecule->setComVel(newComVelOrder4[i]);

            // if(molecule->getMolecularStructureName()=="N2"){
            //     CollisionModel::Atom::rotate2D(angularVelocity*dt, nitrogenOne);
            //     CollisionModel::Atom::rotate2D(angularVelocity*dt, nitrogenTwo);
            //     molecule->getAtoms().at(0)->setRelativePosition(nitrogenOne);
            //     molecule->getAtoms().at(1)->setRelativePosition(nitrogenTwo);
            //     molecule->setAngles(nitrogenAngles);
            // }
            i++;
        }



        steps++;
        dt = dt * globalDelta;

        size_t index = 0;
        for(size_t b = 0; b < nMolecules; ++b){
            for(size_t z = b+1; z < nMolecules; ++z){
                if((moleculesPtr[z]->getComPos() - moleculesPtr[b]->getComPos()).magnitude() > startDistances[index++]){
                    if(trajectoryRecordingActive_ == true && moleculesPtr[z]->getMolecularStructureName() == collisionMolecule_){
                        writeTrajectory((moleculesPtr[z]->getComPos() - moleculesPtr[b]->getComPos()).magnitude(),
                                        moleculesPtr[z]->getComPos(), moleculesPtr[z]->getComVel(), forceMolecules, true, trajectoryOutputStream_.get(), integrationTimeSum, dt);
                    }
                    return wasHit;
                }
                if((moleculesPtr[z]->getComPos() - moleculesPtr[b]->getComPos()).magnitude() <= requiredRad){
                    wasHit=true;
                }
            }
        }
    }
  
    return false;
}
