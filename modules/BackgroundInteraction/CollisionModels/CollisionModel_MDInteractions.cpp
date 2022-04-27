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

CollisionModel::MDInteractionsModel::MDInteractionsModel(double staticPressure,
                                                        double staticTemperature,
                                                        double collisionGasMassAmu,
                                                        double collisionGasDiameterM, 
                                                        double collisionGasPolarizabilityM3) :
        MDInteractionsModel(
        getConstantDoubleFunction(staticPressure),
        getConstantVectorFunction(Core::Vector(0.0, 0.0, 0.0)),
        staticTemperature,
        collisionGasMassAmu,
        collisionGasDiameterM,
        collisionGasPolarizabilityM3) { }

CollisionModel::MDInteractionsModel::MDInteractionsModel(std::function<double(Core::Vector& location)> pressureFunction,
                                                        std::function<Core::Vector(Core::Vector& location)> velocityFunction,
                                                        double staticTemperature,
                                                        double collisionGasMassAmu,
                                                        double collisionGasDiameterM, 
                                                        double collisionGasPolarizabilityM3) :
        MDInteractionsModel(
                std::move(pressureFunction),
                std::move(velocityFunction),
                getConstantDoubleFunction(staticTemperature),
                collisionGasMassAmu,
                collisionGasDiameterM, 
                collisionGasPolarizabilityM3) { }

CollisionModel::MDInteractionsModel::MDInteractionsModel(std::function<double(Core::Vector& location)> pressureFunction,
                                                        std::function<Core::Vector(Core::Vector& location)> velocityFunction,
                                                        std::function<double(const Core::Vector&)> temperatureFunction,
                                                        double collisionGasMassAmu,
                                                        double collisionGasDiameterM, 
                                                        double collisionGasPolarizabilityM3) :
        collisionGasMass_kg_(collisionGasMassAmu*Core::AMU_TO_KG),
        collisionGasDiameter_m_(collisionGasDiameterM),
        collisionGasPolarizability_m3_(collisionGasPolarizabilityM3),
        pressureFunction_(std::move(pressureFunction)),
        velocityFunction_(std::move(velocityFunction)),
        temperatureFunction_(std::move(temperatureFunction)) { }

void CollisionModel::MDInteractionsModel::initializeModelParameters(CollisionModel::Molecule& /*mole*/) const {

}

void CollisionModel::MDInteractionsModel::updateModelParameters(CollisionModel::Molecule& /*mole*/) const {

}

void CollisionModel::MDInteractionsModel::modifyAcceleration(Core::Vector& /*acceleration*/, CollisionModel::Molecule& /*mole*/, double /*dt*/) {
    
}

void CollisionModel::MDInteractionsModel::modifyVelocity(CollisionModel::Molecule& mole, CollisionModel::Molecule& bgMole, double dt) {

    Core::RandomSource* rndSource = Core::globalRandomGeneratorPool->getThreadRandomSource();

    // Calculate collision cross section between particle and collision gas:
    // FIXME: add diamater function for molecule (let it act like a point particle before collision)
    double sigma_m2 = M_PI * std::pow( (mole.getDiameter() + collisionGasDiameter_m_)/2.0, 2.0);
    Core::Vector moleculeComPosition = mole.getComPos();
    double localPressure_Pa = pressureFunction_(moleculeComPosition);

    if (Core::isDoubleEqual(localPressure_Pa, 0.0)){
        return; //pressure 0 means no collision at all
    }

    // Transform the frame of reference in a frame where the mean background gas velocity is zero.
    Core::Vector vGasMean = velocityFunction_(moleculeComPosition);
    Core::Vector vFrameMeanBackRest = mole.getComVel() - vGasMean;

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

    // Decide if a collision actually happens:
    if (rndSource->uniformRealRndValue() > collisionProb){
        return; // no collision takes place
    }           

    // Collision happens
    // TODO: Construct the actual molecule and its atoms - do this before call to this function 
    // TODO: Construct the background gas particle - do this before call to this function

    // Switch to COM frame 

    Core::Vector momentumSum = Core::Vector(0.0, 0.0, 0.0);
    Core::Vector positionSum = Core::Vector(0.0, 0.0, 0.0);
    double massSum = 0;
    std::vector<CollisionModel::Molecule*> moleculesPtr = {&mole, &bgMole};
    for(auto* molecule : moleculesPtr){
        momentumSum += molecule->getComVel() * molecule->getMass();
        positionSum += molecule->getComPos() * molecule->getMass();
        massSum += molecule->getMass();
    }
    for(auto* molecule : moleculesPtr){
        molecule->setComVel(molecule->getComVel() - (momentumSum / massSum));
        molecule->setComPos(molecule->getComPos() - (positionSum / massSum));
    }

    // Call the sub-integrator
    double finalTime = 4E-13; //  final integration time in seconds
    double timeStep = 1E-15; // step size in seconds 
    leapfrogIntern(moleculesPtr, timeStep, finalTime);


}

void CollisionModel::MDInteractionsModel::modifyPosition(Core::Vector& /*position*/, CollisionModel::Molecule& /*mole*/, double /*dt*/) {

}

void CollisionModel::MDInteractionsModel::leapfrogIntern(std::vector<CollisionModel::Molecule*> moleculesPtr, double dt, double finalTime){

    int nSteps = int(round(finalTime/dt));
    //size_t nMolecules = moleculesPtr.size();

    std::vector<Core::Vector> forceMolecules = forceFieldMD(moleculesPtr);

    // do the first half step for the velocity, as per leapfrog definition
    size_t i = 0;
    for(auto* molecule : moleculesPtr){
        Core::Vector newComVel =  molecule->getComVel() + forceMolecules.at(i) / molecule->getMass() * dt/2;
        molecule->setComVel(newComVel);
        i++;
    }

    // start the actual leapfrog iteration
    for (int j = 0; j < nSteps; j++){
        
        // time step for the new position 
        for(auto* molecule : moleculesPtr){
            Core::Vector newComPos =  molecule->getComPos() + molecule->getComVel() * dt;
            molecule->setComPos(newComPos);
        }

        // recalculate the force
        std::vector<Core::Vector> forceMolecules = forceFieldMD(moleculesPtr);  
        i = 0;
        // time step for the new velocity
        for(auto* molecule : moleculesPtr){
            Core::Vector newComVel =  molecule->getComVel() + forceMolecules.at(i) / molecule->getMass() * dt;
            molecule->setComVel(newComVel);
            i++;
        }
    }

}

std::vector<Core::Vector> CollisionModel::MDInteractionsModel::forceFieldMD(std::vector<CollisionModel::Molecule*> moleculesPtr){

    size_t nMolecules = moleculesPtr.size();
    std::vector<Core::Vector> forceMolecules; // save all the forces acting on each molecule

    // each molecule interacts with each other molecule
    for (size_t i = 0; i < nMolecules; i++){
        for (size_t j = i+1; j < nMolecules; j++){
            /* 
            * therefore we need the interaction between each atom of a molecule with the atoms of the
            * other one 
            */ 

            // construct E-field acting on the molecule
            std::array<double, 3> eField;
            std::array<double, 6> eFieldDerivative;

            for(auto* atomI : moleculesPtr.at(i)->getAtoms()){
                for(auto* atomJ : moleculesPtr.at(j)->getAtoms()){
                    
                    // First contribution: Lennard-Jones potential 
                    // This always contributes to the experienced force 
                    Core::Vector absPosAtomI = moleculesPtr.at(i)->getComPos() + atomI->getRelativePosition();
                    Core::Vector absPosAtomJ = moleculesPtr.at(j)->getComPos() + atomJ->getRelativePosition();

                    Core::Vector distance = absPosAtomI - absPosAtomJ;
                    double distanceSquared = distance.magnitudeSquared();
                    double distanceSquaredInverse = 1/distanceSquared;
                    double sigma6 = pow(CollisionModel::Atom::calcLJSig(*atomI, *atomJ), 6);
                    double epsilon = CollisionModel::Atom::calcLJEps(*atomI, *atomJ);
                    double ljFactor = 24 * epsilon * pow(distanceSquaredInverse, 4) * 
                                        (2 * pow(distanceSquaredInverse, 3) * sigma6 * sigma6 - sigma6);

                    // calculate the force that acts on the atoms and add it to the overall force on the molecule
                    Core::Vector atomForce;
                    atomForce.x(distance.x() * ljFactor);
                    atomForce.y(distance.y() * ljFactor);
                    atomForce.z(distance.z() * ljFactor);
                    forceMolecules.at(i) += atomForce;
                    forceMolecules.at(j) += atomForce * (-1); 

                    // Second contribution: C4 ion-induced dipole potential
                    // This requires an ion and one neutrally charged molecule to be present
                    double distanceCubed = distanceSquared * sqrt(distanceSquared);
                    double currentCharge = 0;
                    // Check if one of the molecules is an ion and the other one is not
                    if(moleculesPtr.at(i)->getIsIon() == true && 
                        moleculesPtr.at(j)->getIsIon() == false){
                        currentCharge = atomI->getCharge();

                    }else if (moleculesPtr.at(i)->getIsIon() == false && 
                                moleculesPtr.at(j)->getIsIon() == true){
                        currentCharge = atomJ->getCharge();
                    }

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

                    // Third contribution: ion <-> permanent dipole potential
                    // This requires an ion and a dipole to be present 
                    if(moleculesPtr.at(i)->getIsIon() == true && 
                        moleculesPtr.at(j)->getIsDipole() == true){

                        currentCharge = atomI->getCharge();
                        double dipoleDistanceScalar =  moleculesPtr.at(j)->getDipole().x() * distance.x() + 
                                                        moleculesPtr.at(j)->getDipole().y() * distance.y() + 
                                                        moleculesPtr.at(j)->getDipole().z() * distance.z();
                        Core::Vector ionDipoleForce;
                        ionDipoleForce.x(-currentCharge * 1/Core::ELECTRIC_CONSTANT * 
                                            (1/distanceCubed * moleculesPtr.at(j)->getDipole().x() - 
                                            3 * dipoleDistanceScalar * 1/(distanceCubed*distanceSquared) * distance.x()) );
                        ionDipoleForce.y(-currentCharge * 1/Core::ELECTRIC_CONSTANT * 
                                            (1/distanceCubed * moleculesPtr.at(j)->getDipole().y() - 
                                            3 * dipoleDistanceScalar * 1/(distanceCubed*distanceSquared) * distance.y()) );
                        ionDipoleForce.z(-currentCharge * 1/Core::ELECTRIC_CONSTANT * 
                                            (1/distanceCubed * moleculesPtr.at(j)->getDipole().z() - 
                                            3 * dipoleDistanceScalar * 1/(distanceCubed*distanceSquared) * distance.z()) );
                        forceMolecules.at(i) += ionDipoleForce;
                        forceMolecules.at(j) += ionDipoleForce * (-1); 

                    }else if (moleculesPtr.at(i)->getIsDipole() == true && 
                                moleculesPtr.at(j)->getIsIon() == true){

                        currentCharge = atomJ->getCharge();
                        double dipoleDistanceScalar =  moleculesPtr.at(i)->getDipole().x() * distance.x() + 
                                                        moleculesPtr.at(i)->getDipole().y() * distance.y() + 
                                                        moleculesPtr.at(i)->getDipole().z() * distance.z();
                        Core::Vector ionDipoleForce;
                        ionDipoleForce.x(-currentCharge * 1/Core::ELECTRIC_CONSTANT * 
                                            (1/distanceCubed * moleculesPtr.at(i)->getDipole().x() - 
                                            3 * dipoleDistanceScalar * 1/(distanceCubed*distanceSquared) * distance.x()) );
                        ionDipoleForce.y(-currentCharge * 1/Core::ELECTRIC_CONSTANT * 
                                            (1/distanceCubed * moleculesPtr.at(i)->getDipole().y() - 
                                            3 * dipoleDistanceScalar * 1/(distanceCubed*distanceSquared) * distance.y()) );
                        ionDipoleForce.z(-currentCharge * 1/Core::ELECTRIC_CONSTANT * 
                                            (1/distanceCubed * moleculesPtr.at(i)->getDipole().z() - 
                                            3 * dipoleDistanceScalar * 1/(distanceCubed*distanceSquared) * distance.z()) );
                        forceMolecules.at(i) += ionDipoleForce;
                        forceMolecules.at(j) += ionDipoleForce * (-1); 
                    }
                }
            }
            // add the C4 ion induced force
            Core::Vector ionInducedForce;
            ionInducedForce.x(1/Core::ELECTRIC_CONSTANT * collisionGasPolarizability_m3_ * 
                                (eField[0]*eFieldDerivative[0] + eField[1]*eFieldDerivative[1] + eField[2]*eFieldDerivative[5]));
            ionInducedForce.y(1/Core::ELECTRIC_CONSTANT * collisionGasPolarizability_m3_ * 
                                (eField[0]*eFieldDerivative[1] + eField[1]*eFieldDerivative[2] + eField[2]*eFieldDerivative[3]));
            ionInducedForce.z(1/Core::ELECTRIC_CONSTANT * collisionGasPolarizability_m3_ * 
                                (eField[0]*eFieldDerivative[5] + eField[1]*eFieldDerivative[3] + eField[2]*eFieldDerivative[4]));
            forceMolecules.at(i) += ionInducedForce;
            forceMolecules.at(j) += ionInducedForce * (-1);
        }
    }
    return forceMolecules;
}

void CollisionModel::MDInteractionsModel::initializeModelParameters(Core::Particle& /*ion*/) const {

}

void CollisionModel::MDInteractionsModel::updateModelParameters(Core::Particle& /*ion*/) const {

}

void CollisionModel::MDInteractionsModel::modifyAcceleration(Core::Vector& /*acceleration*/, Core::Particle& /*particle*/,
                                                         double /*dt*/) {

}

void CollisionModel::MDInteractionsModel::modifyVelocity(Core::Particle& /*particle*/, double /*dt*/) {

}

void CollisionModel::MDInteractionsModel::modifyPosition(Core::Vector& /*position*/, Core::Particle& /*particle*/, double /*dt*/) {

}