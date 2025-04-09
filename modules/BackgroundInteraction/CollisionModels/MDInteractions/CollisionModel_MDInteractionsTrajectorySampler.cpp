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

#include "CollisionModel_MDInteractionsTrajectorySampler.hpp"
#include "Core_math.hpp"
#include "Core_randomGenerators.hpp"
#include <cmath>
#include <array>
#include <fstream>

CollisionModel::MDInteractionsTrajectorySampler::MDInteractionsTrajectorySampler(
            std::string collisionMolecule,
            double integrationTime,
            double subTimeStep,
            std::unique_ptr<CollisionModel::AbstractMDForceField> forceField,
            std::unordered_map<std::string, std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection,
            Core::Vector startPosition,
            Core::Vector startVelocity) :
        collisionMolecule_(collisionMolecule),
        integrationTime_(integrationTime),
        subTimeStep_(subTimeStep),
        forceField_(std::move(forceField)),
        molecularStructureCollection_(std::move(molecularStructureCollection)),
        startPosition_(startPosition), 
        startVelocity_(startVelocity){}
/**
 * Activates trajectory writing and sets trajectory writer configuration
 * @param trajectoryFileName
 * @param trajectoryDistance
 */
void CollisionModel::MDInteractionsTrajectorySampler::setTrajectoryWriter(const std::string& trajectoryFileName,
                                                              double minimalSampleInterval) {

    trajectoryWriter_= std::make_unique<CollisionModel::MDTrajectoryWriter>(trajectoryFileName, minimalSampleInterval);
}


void CollisionModel::MDInteractionsTrajectorySampler::calculateTrajectory(Core::Particle& particle) {

    int iterations = 0;

    // Construct the actual molecule and its atoms
    Molecule mole(Core::Vector(0.0, 0.0, 0.0), Core::Vector(0.0, 0.0, 0.0), particle.getMolecularStructure());

    // Construct the background gas particle
    Molecule bgMole(startPosition_, startVelocity_,
                                        molecularStructureCollection_.at(collisionMolecule_));

    // TODO: Rotate background gas molecule here
    //moleculesPtr[1]->setAngles(nitrogenAngles);
    std::vector<CollisionModel::Molecule*> moleculesPtr = {&mole, &bgMole};

    // possible check for energy conservation
    std::vector<Core::Vector> startVelocity;
    double startEnergy = 0;
    for(auto* molecule : moleculesPtr){
        startVelocity.push_back(molecule->getComVel());
        startEnergy += 0.5 * molecule->getMass() * molecule->getComVel().magnitudeSquared();
    }

    // Call the sub-integrator
    double finalTime = integrationTime_; // final integration time in seconds
    double timeStep = subTimeStep_; // step size in seconds


    //trajectorySuccess = rk4Intern(moleculesPtr, timeStep, finalTime, collisionRadius);
    //trajectorySuccess = leapfrogIntern(moleculesPtr, timeStep, finalTime, collisionRadius);
    bool trajectorySuccess = rk4InternAdaptiveStep(moleculesPtr, timeStep, finalTime);

    double endEnergy = 0;
    for(auto* molecule : moleculesPtr){
        endEnergy += 0.5 * molecule->getMass() * molecule->getComVel().magnitudeSquared();
    }
    particle.setVelocity(mole.getComVel() + particle.getVelocity());
    ++iterations;
}

bool CollisionModel::MDInteractionsTrajectorySampler::leapfrogIntern(std::vector<CollisionModel::Molecule*> moleculesPtr, double dt, double finalTime, double requiredRad){
    bool wasHit = false;
    double distance = 0.0;
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
        //if(molecule->getMolecularStructureName() == collisionMolecule_){
            molecule->setComVel(newComVel);
        //}
        i++;
    }

    // start the actual leapfrog iteration
    for (int j = 0; j < nSteps; j++){
        // time step for the new position
        i = 0;
        double energyEnd = 0;

        for(size_t k = 0; k < moleculesPtr_size; ++k){
            for(size_t l = k+1; l < moleculesPtr_size; ++l){
                distance = (moleculesPtr[l]->getComPos() - moleculesPtr[k]->getComPos()).magnitude();
            }
        }

        for(auto* molecule : moleculesPtr){
            if(trajectoryWriter_ != nullptr && molecule->getMolecularStructureName() == collisionMolecule_){
                trajectoryWriter_->writeTrajectorySample(
                    j * dt, dt, molecule->getComPos(), molecule->getComVel(),
                    moleculesPtr[0]->getComPos(),forceMolecules, distance);
            }
            Core::Vector newComPos =  molecule->getComPos() + molecule->getComVel() * dt;
            //if(molecule->getMolecularStructureName() == collisionMolecule_){
                molecule->setComPos(newComPos);
            //}
            energyEnd += 0.5 * molecule->getComVel().magnitudeSquared() * molecule->getMass();
            i++;
        }

        size_t index = 0;
        for(size_t b = 0; b < moleculesPtr_size; ++b){
            for(size_t z = b+1; z < moleculesPtr_size; ++z){
                if((moleculesPtr.at(z)->getComPos() - moleculesPtr.at(b)->getComPos()).magnitude() > startDistances.at(index++)){
                    if(trajectoryWriter_ != nullptr && moleculesPtr[z]->getMolecularStructureName() == collisionMolecule_ && (j+1) == nSteps){
                        double distance = (moleculesPtr[z]->getComPos() - moleculesPtr[b]->getComPos()).magnitude();
                        trajectoryWriter_->writeTrajectorySample(
                            j*dt, dt, moleculesPtr[z]->getComPos(), moleculesPtr[z]->getComVel(),
                            moleculesPtr[0]->getComPos(),forceMolecules, distance);
                    }
                    //return wasHit;
                }
                else if((moleculesPtr.at(z)->getComPos() - moleculesPtr.at(b)->getComPos()).magnitude() <= requiredRad){
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
            //if(molecule->getMolecularStructureName() == collisionMolecule_){
                molecule->setComVel(newComVel);
            //}
            i++;
        }
    }
    return false;
}

bool CollisionModel::MDInteractionsTrajectorySampler::rk4Intern(std::vector<CollisionModel::Molecule*> moleculesPtr, double dt, double finalTime,
                                                                    double requiredRad){

    int nSteps = int(round(finalTime/dt));
    size_t nMolecules = moleculesPtr.size();
    std::vector<Core::Vector> forceMolecules(nMolecules);

    bool wasHit = false;
    double distance = 0.0;
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
            if(molecule->getMolecularStructureName() == collisionMolecule_){
                velocityMolecules.at(i) = molecule->getComVel();
                positionMolecules.at(i) = molecule->getComPos();
            }else{
                velocityMolecules.at(i) = {0,0,0};
                positionMolecules.at(i) = {0,0,0};
            }
            i++;
        }
        
        std::vector<Core::Vector> initialPositionMolecules(nMolecules);
        std::vector<Core::Vector> initialVelocityMolecules(nMolecules);
        for(size_t k = 0; k < nMolecules; k++){
            initialPositionMolecules.at(k) = Core::Vector( positionMolecules.at(k).x(), positionMolecules.at(k).y(),positionMolecules.at(k).z() );
            initialVelocityMolecules.at(k) = Core::Vector( velocityMolecules.at(k).x(), velocityMolecules.at(k).y(), velocityMolecules.at(k).z() );
        }
        
        for(size_t k = 0; k < nMolecules; ++k){
            for(size_t l = k+1; l < nMolecules; ++l){
                distance = (moleculesPtr[l]->getComPos() - moleculesPtr[k]->getComPos()).magnitude();
            }
        }

        

        double length[3] = {1./2, 1./2, 1};
        std::vector<double> mass;
        i = 0;
        for(auto* molecule : moleculesPtr){
            mass.push_back(molecule->getMass());
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
            Core::Vector newComVel = initialVelocityMolecules.at(i) + (k[0][i]+ k[1][i]*2 + k[2][i]*2 + k[3][i]) * 1./6;
            if(trajectoryWriter_ != nullptr && molecule->getMolecularStructureName() == collisionMolecule_){
                trajectoryWriter_->writeTrajectorySample(
                    j*dt, dt, molecule->getComPos(), molecule->getComVel(),
                    moleculesPtr[0]->getComPos(),forceMolecules, distance);
                }
            if(molecule->getMolecularStructureName() == collisionMolecule_){
                molecule->setComPos(newComPos);
                molecule->setComVel(newComVel);
            }
            i++;
        }

        size_t index = 0;
        for(size_t k = 0; k < nMolecules; ++k){
            for(size_t l = k+1; l < nMolecules; ++l){
                if((moleculesPtr[l]->getComPos() - moleculesPtr[k]->getComPos()).magnitude() > startDistances[index++]){
                    if(trajectoryWriter_ != nullptr && moleculesPtr[l]->getMolecularStructureName() == collisionMolecule_ && (j+1) == nSteps){
                        double distance = (moleculesPtr[l]->getComPos() - moleculesPtr[k]->getComPos()).magnitude();
                        trajectoryWriter_->writeTrajectorySample(
                            j*dt, dt, moleculesPtr[l]->getComPos(), moleculesPtr[l]->getComVel(),
                            moleculesPtr[0]->getComPos(),forceMolecules, distance);
                    }
                    //return wasHit;
                }
                if((moleculesPtr[l]->getComPos() - moleculesPtr[k]->getComPos()).magnitude() <= requiredRad){
                    wasHit=true;
                }
            }
        }
    }
    return false;
}


bool CollisionModel::MDInteractionsTrajectorySampler::rk4InternAdaptiveStep(std::vector<CollisionModel::Molecule*> moleculesPtr, double dt, double finalTime){

    double integrationTimeSum = 0;
    size_t nMolecules = moleculesPtr.size();
    std::vector<Core::Vector> forceMolecules(nMolecules);

    size_t i = 0;
    int steps = 0;
    double distance = 0.0;

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
            initialPositionMolecules[k] = Core::Vector( positionMolecules[k].x(), positionMolecules[k].y(), positionMolecules[k].z() );
            initialVelocityMolecules[k] = Core::Vector( velocityMolecules[k].x(), velocityMolecules[k].y(), velocityMolecules[k].z() );
        }

        i = 0;
        for(auto* molecule : moleculesPtr){
            mass[i] = molecule->getMass();
            i++;
        }
        
        forceField_->calculateForceField(moleculesPtr, forceMolecules);

        for(size_t q = 0; q < nMolecules; q++){
            k[0][q] = forceMolecules[q] * dt / mass[q];
            l[0][q] = velocityMolecules[q] * dt;
        }

        for(size_t n = 1; n < 6; n++){
            i = 0;
            for(size_t i=0; i<moleculesPtr.size(); ++i){
                positionMolecules[i] = initialPositionMolecules[i];
            }
            for(size_t m = 0; m < 6; m++){
                i = 0;
                for(auto* molecule : moleculesPtr){
                    positionMolecules[i] += l[m][i]*weight[n-1][m];
                    //TODO: Add Switch statement for ion movement:|| integrate all = true
                    if(molecule->getMolecularStructureName() == collisionMolecule_){
                        molecule->setComPos(positionMolecules[i]);
                    }
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
        std::array<double,2> R;

        for(size_t p = 0; p < 2; p++){
            if(fabs(newComVelOrder5[p].magnitude()) != 0)
                R[p] = fabs(newComVelOrder4[p].magnitude()-newComVelOrder5[p].magnitude())/fabs(newComVelOrder5[p].magnitude());
            else
                R[p] = 0;
        }

        double globalR = std::max({R[0],R[1]});
        double tolerance = 1e-8;

        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wfloat-equal"
        if (globalR == 0){
            globalR = 1e-15;
        }
        #pragma GCC diagnostic pop
        double globalDelta = 0.84 * std::pow((tolerance/globalR), 1./4);
        
        i = 0;
        integrationTimeSum += dt;
        for(auto* molecule : moleculesPtr){

            //TODO: Add Switch statement for ion movement:|| integrate all = true
            if(molecule->getMolecularStructureName() == collisionMolecule_){
                molecule->setComPos(newComPosOrder4[i]);
                molecule->setComVel(newComVelOrder4[i]); 
            }
            if(trajectoryWriter_ != nullptr && molecule->getMolecularStructureName() == collisionMolecule_ && integrationTimeSum < finalTime){
                trajectoryWriter_->writeTrajectorySample(
                    integrationTimeSum, dt, molecule->getComPos(), molecule->getComVel(),
                    moleculesPtr[0]->getComPos(),forceMolecules, distance);
            }
            i++;
        }
        steps++;
        dt = dt * globalDelta;
    }
    if(trajectoryWriter_ != nullptr) {
        trajectoryWriter_->writeTrajectoryDelimiter();
    }
    return true;
}

