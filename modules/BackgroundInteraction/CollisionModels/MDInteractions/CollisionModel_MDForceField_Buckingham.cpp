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
 ****************************/

#include "CollisionModel_MDForceField_Buckingham.hpp"
#include "Core_math.hpp"
#include <array>
#include <iostream>
#include <cmath>

CollisionModel::MDForceField_Buckingham::MDForceField_Buckingham(double collisionGasPolarizability_m3, 
                                                                    std::string potentials, 
                                                                    bool rotActive):
    collisionGasPolarizability_m3_(collisionGasPolarizability_m3), 
    potentialsFF_(potentials), 
    rotActive_(rotActive)
{}

void CollisionModel::MDForceField_Buckingham::calculateForceField(std::vector<CollisionModel::Molecule*>& moleculesPtr,
                                                                std::vector<Core::Vector>& forceMolecules, 
                                                                std::vector<Core::Vector>& torqueMolecules) {

    // save all the forces acting on each molecule
    CollisionModel::Molecule* ion = moleculesPtr[0];
    CollisionModel::Molecule* bgGas = moleculesPtr[1];
    forceMolecules[0] = Core::Vector(0.0, 0.0, 0.0);
    forceMolecules[1] = Core::Vector(0.0, 0.0, 0.0);
    bool isN2 = false;
    bool isN2Approx = false;
    bool isCO2 = false;

    if(bgGas->getMolecularStructureName()=="N2"){
        isN2 = true;
    }else if(bgGas->getMolecularStructureName()=="N2Approx"){
        isN2Approx = true;
    }else if(bgGas->getMolecularStructureName()=="CO2"){
        isCO2 = true;
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
            Core::Vector absPosAtomI = ion->genWorldFramePosition(atomI);
            Core::Vector absPosAtomJ = bgGas->genWorldFramePosition(atomJ);
            Core::Vector relPosAtomI = ion->genWorldFramePositionRelativeCOM(atomI);
            Core::Vector relPosAtomJ = bgGas->genWorldFramePositionRelativeCOM(atomJ);

            Core::Vector distance = absPosAtomI - absPosAtomJ;
            double distanceSquared = distance.magnitudeSquared();
            double distanceAbs = sqrt(distanceSquared);
            double distanceSquaredInverse = 1./distanceSquared;
            double distanceCubed = distanceSquared * sqrt(distanceSquared);
            double rmax = interactionMap.get(*atomI, *atomJ);

            if(distance.magnitude() > 100e-10){
                return;
            }
            
            if(potentialsFF_ == "ALL" || potentialsFF_ == "VDW" || potentialsFF_ == "VDWII" || potentialsFF_ == "VDWQUAD"){
                double sigma = CollisionModel::Atom::calcLJSig(*atomI, *atomJ);
                double sigma6 = sigma * sigma * sigma * sigma * sigma * sigma;
                double epsilon = CollisionModel::Atom::calcLJEps(*atomI, *atomJ);
                double ljFactor = 0;
                if(distanceAbs <= rmax){
                    ljFactor = (-10) * epsilon * 1/distanceAbs * -1.84e5 * 12 * exp(-12*distanceAbs/sigma) * 1/sigma;
                }else{
                    ljFactor = (-1) * epsilon * 1/distanceAbs * (-1.84e5 * 12 * exp(-12*distanceAbs/sigma) * 1/sigma + 
                                            2.25 * 6 * sigma6 * distanceSquaredInverse*distanceSquaredInverse*distanceSquaredInverse * 1/distanceAbs);
                }
                    
                // calculate the force that acts on the atoms and add it to the overall force on the molecule
                Core::Vector atomForce;
                atomForce.x(distance.x() * ljFactor);
                atomForce.y(distance.y() * ljFactor);
                atomForce.z(distance.z() * ljFactor);

                forceMolecules[0] += atomForce;
                forceMolecules[1] += atomForce * (-1);
                 if(rotActive_){

                    std::vector<Core::Vector> positions = {relPosAtomI, relPosAtomJ};
                    calculateTorque(positions, atomForce, torqueMolecules);
                }
                
            }

            // Second contribution: C4 ion-induced dipole potential
            // This requires an ion and one neutrally charged molecule to be present
            if(potentialsFF_ == "ALL" || potentialsFF_ == "II" || potentialsFF_ == "VDWII" || potentialsFF_ == "VDWIIQUAD"){
                
                double currentCharge = 0;
                // Check if one of the molecules is an ion and the other one is not
                if(isN2 || isCO2){
                    if(int(ceil(fabs(atomI->getCharge()/Core::ELEMENTARY_CHARGE))) != 0 &&
                            atomJ->getType() == CollisionModel::Atom::AtomType::COM){

                        currentCharge = atomI->getCharge();
                    }

                }else{
                    if(int(ceil(fabs(atomI->getCharge()/Core::ELEMENTARY_CHARGE))) != 0 &&
                            atomJ->getType() != CollisionModel::Atom::AtomType::COM){

                        currentCharge = atomI->getCharge();

                    }
                }

                eField[0] = distance.x() * currentCharge / distanceCubed; // E-field in x
                eField[1] = distance.y() * currentCharge / distanceCubed; // E-field in y
                eField[2] = distance.z() * currentCharge / distanceCubed; // E-field in z

                // derivative x to x
                eFieldDerivative[0] = currentCharge / distanceCubed -
                        3 * currentCharge * distance.x() * distance.x() / (distanceCubed * distanceSquared);
                // derivative x to y
                eFieldDerivative[1] = -3 * currentCharge * distance.x() * distance.y() / (distanceCubed * distanceSquared);
                // derivative y to y
                eFieldDerivative[2] = currentCharge / distanceCubed -
                        3 * currentCharge * distance.y() * distance.y() / (distanceCubed * distanceSquared);
                // derivative y to z
                eFieldDerivative[3] = -3 * currentCharge * distance.y() * distance.z() / (distanceCubed * distanceSquared);
                // derivative z to z
                eFieldDerivative[4] = currentCharge / distanceCubed -
                        3 * currentCharge * distance.z() * distance.z() / (distanceCubed * distanceSquared);
                // derivative x to z
                eFieldDerivative[5] = -3 * currentCharge * distance.x() * distance.z() / (distanceCubed * distanceSquared);

                Core::Vector ionInducedForce;
                if(isN2Approx){
                    double xIonInducedForce = 1./(Core::ELECTRIC_CONSTANT) * collisionGasPolarizability_m3_/2 *
                            (eField[0]*eFieldDerivative[0] + eField[1]*eFieldDerivative[1] + eField[2]*eFieldDerivative[5]);
                    double yIonInducedForce = 1./(Core::ELECTRIC_CONSTANT) * collisionGasPolarizability_m3_/2 *
                            (eField[0]*eFieldDerivative[1] + eField[1]*eFieldDerivative[2] + eField[2]*eFieldDerivative[3]);
                    double zIonInducedForce = 1./(Core::ELECTRIC_CONSTANT) * collisionGasPolarizability_m3_/2 *
                            (eField[0]*eFieldDerivative[5] + eField[1]*eFieldDerivative[3] + eField[2]*eFieldDerivative[4]);
                    ionInducedForce.x(xIonInducedForce);
                    ionInducedForce.y(yIonInducedForce);
                    ionInducedForce.z(zIonInducedForce);
                }else{
                    double xIonInducedForce = 1./(Core::ELECTRIC_CONSTANT) * collisionGasPolarizability_m3_ *
                            (eField[0]*eFieldDerivative[0] + eField[1]*eFieldDerivative[1] + eField[2]*eFieldDerivative[5]);
                    double yIonInducedForce = 1./(Core::ELECTRIC_CONSTANT) * collisionGasPolarizability_m3_ *
                            (eField[0]*eFieldDerivative[1] + eField[1]*eFieldDerivative[2] + eField[2]*eFieldDerivative[3]);
                    double zIonInducedForce = 1./(Core::ELECTRIC_CONSTANT) * collisionGasPolarizability_m3_ *
                            (eField[0]*eFieldDerivative[5] + eField[1]*eFieldDerivative[3] + eField[2]*eFieldDerivative[4]);
                    ionInducedForce.x(xIonInducedForce);
                    ionInducedForce.y(yIonInducedForce);
                    ionInducedForce.z(zIonInducedForce);
                }
                double cutoff = 1.5*CollisionModel::Atom::calcLJSig(*atomI, *atomJ);
                if(distance.magnitude() < cutoff && rotActive_){
                    
                    double dist = distance.magnitude();  
                    double a = ((dist-cutoff)*1e10);
                    double decayFuncValue = (1./(1+exp(a*a)*exp(a*a))*2);
                    ionInducedForce = ionInducedForce * decayFuncValue; 
                }

                forceMolecules[0] += ionInducedForce;
                forceMolecules[1] += ionInducedForce * (-1);

                if(rotActive_){

                    std::vector<Core::Vector> positions = {relPosAtomI, relPosAtomJ};
                    calculateTorque(positions, ionInducedForce, torqueMolecules);
                }
            }
            



            // Third contribution: ion <-> permanent dipole potential
            // This requires an ion and a dipole to be present
            // double dipoleDistanceScalar = 0;
            //double dipoleX = 0, dipoleY = 0, dipoleZ = 0;
            if(potentialsFF_ == "ALL" || potentialsFF_ == "PDIPOLE" || potentialsFF_ == "VDWPDIPOLE"){

                double currentCharge = 0;
                if(int(atomI->getCharge()/Core::ELEMENTARY_CHARGE) != 0 &&
                        moleculesPtr[1]->getIsDipole() == true){

                    currentCharge = atomI->getCharge();
                    /*dipoleX = moleculesPtr[1]->getDipole().x();
                    dipoleY = moleculesPtr[1]->getDipole().y();
                    dipoleZ = moleculesPtr[1]->getDipole().z();
                    dipoleDistanceScalar =  dipoleX * distance.x() +
                            dipoleY * distance.y() +
                            dipoleZ * distance.z();*/

                }else if (moleculesPtr[0]->getIsDipole() == true &&
                        int(atomJ->getCharge()/Core::ELEMENTARY_CHARGE) != 0){

                    currentCharge = atomJ->getCharge();
                    /*dipoleX = moleculesPtr[0]->getDipole().x();
                    dipoleY = moleculesPtr[0]->getDipole().y();
                    dipoleZ = moleculesPtr[0]->getDipole().z();
                    dipoleDistanceScalar =  dipoleX * distance.x() +
                            dipoleY * distance.y() +
                            dipoleZ * distance.z();*/
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
            }
            

            // Fourth contribution: quadrupole moment if background gas is N2
            // This requires an ion and N2 to be present
            if(potentialsFF_ == "ALL" || potentialsFF_ == "QUAD" || potentialsFF_ == "VDWQUAD" || potentialsFF_ == "VDWIIQUAD"){
                double partialChargeN2 = 0;
                double currentCharge = 0;
                if(int(ceil(fabs(atomI->getCharge()/Core::ELEMENTARY_CHARGE))) != 0 &&
                        (isN2 == true || isN2Approx == true)){

                    currentCharge = atomI->getCharge();
                    partialChargeN2 = atomJ->getPartCharge();

                }else if ((isN2 == true || isN2Approx == true) &&
                        int(ceil(fabs(atomJ->getCharge()/Core::ELEMENTARY_CHARGE))) != 0){

                    currentCharge = atomJ->getCharge();
                    partialChargeN2 = atomI->getPartCharge();
                }
                Core::Vector quadrupoleForce;
                quadrupoleForce.x(currentCharge * partialChargeN2 * 1./Core::ELECTRIC_CONSTANT * distance.x() / distanceCubed );
                quadrupoleForce.y(currentCharge * partialChargeN2 * 1./Core::ELECTRIC_CONSTANT * distance.y() / distanceCubed );
                quadrupoleForce.z(currentCharge * partialChargeN2 * 1./Core::ELECTRIC_CONSTANT * distance.z() / distanceCubed );
                forceMolecules[0] += quadrupoleForce;
                forceMolecules[1] += quadrupoleForce * (-1);
                if(rotActive_){
                    
                    std::vector<Core::Vector> positions = {relPosAtomI, relPosAtomJ};
                    calculateTorque(positions, quadrupoleForce, torqueMolecules);
                }
            }
        }
    }
}

void CollisionModel::MDForceField_Buckingham::calculateForceFieldComponents(std::vector<CollisionModel::Molecule*>& moleculesPtr,
                                                                std::vector<Core::Vector>& forceMolecules, 
                                                                Core::Vector& forceVDW, 
                                                                Core::Vector& forceII) {

    // save all the forces acting on each molecule
    CollisionModel::Molecule* ion = moleculesPtr[0];
    CollisionModel::Molecule* bgGas = moleculesPtr[1];
    forceMolecules[0] = Core::Vector(0.0, 0.0, 0.0);
    forceMolecules[1] = Core::Vector(0.0, 0.0, 0.0);
    forceVDW = Core::Vector(0.0, 0.0, 0.0);
    forceII = Core::Vector(0.0, 0.0, 0.0);
    bool isN2 = false;
    bool isN2Approx = false;
    bool isCO2 = false;

    if(bgGas->getMolecularStructureName()=="N2"){
        isN2 = true;
    }else if(bgGas->getMolecularStructureName()=="N2Approx"){
        isN2Approx = true;
    }else if(bgGas->getMolecularStructureName()=="CO2"){
        isCO2 = true;
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
            Core::Vector absPosAtomI = ion->genWorldFramePosition(atomI);
            Core::Vector absPosAtomJ = bgGas->genWorldFramePosition(atomJ);


            Core::Vector distance = absPosAtomI - absPosAtomJ;
            double distanceSquared = distance.magnitudeSquared();
            double distanceAbs = sqrt(distanceSquared);
            double distanceSquaredInverse = 1./distanceSquared;
            double distanceCubed = distanceSquared * sqrt(distanceSquared);

            if(distance.magnitude() > 100e-10){
                return;
            }
            
            if(potentialsFF_ == "ALL" || potentialsFF_ == "VDW" || potentialsFF_ == "VDWII"){
                double sigma = CollisionModel::Atom::calcLJSig(*atomI, *atomJ);
                double sigma6 = sigma * sigma * sigma * sigma * sigma * sigma;
                double epsilon = CollisionModel::Atom::calcLJEps(*atomI, *atomJ);
                std::cout << sigma <<  " " <<epsilon << std::endl;
                double ljFactor = (-1) * epsilon * 1/distanceAbs * (-1.84e5 * 12 * exp(-12*distanceAbs/sigma) * 1/sigma + 
                                            2.25 * 6 * sigma6 * distanceSquaredInverse*distanceSquaredInverse*distanceSquaredInverse * 1/distanceAbs);

                // calculate the force that acts on the atoms and add it to the overall force on the molecule
                Core::Vector atomForce;
                atomForce.x(distance.x() * ljFactor);
                atomForce.y(distance.y() * ljFactor);
                atomForce.z(distance.z() * ljFactor);
                forceMolecules[0] += atomForce;
                forceMolecules[1] += atomForce * (-1);
                forceVDW = atomForce * (-1);
            }

            // Second contribution: C4 ion-induced dipole potential
            // This requires an ion and one neutrally charged molecule to be present
            if(potentialsFF_ == "ALL" || potentialsFF_ == "II" || potentialsFF_ == "VDWII"){
                double currentCharge = 0;
                // Check if one of the molecules is an ion and the other one is not
                if(isN2 || isCO2){
                    if(int(ceil(fabs(atomI->getCharge()/Core::ELEMENTARY_CHARGE))) != 0 &&
                            atomJ->getType() == CollisionModel::Atom::AtomType::COM){

                        currentCharge = atomI->getCharge();

                    }else if (int(ceil(fabs(atomJ->getCharge()/Core::ELEMENTARY_CHARGE))) != 0 &&
                            atomI->getType() == CollisionModel::Atom::AtomType::COM){

                        currentCharge = atomJ->getCharge();

                    }
                }else{
                    if(int(ceil(fabs(atomI->getCharge()/Core::ELEMENTARY_CHARGE))) != 0 &&
                            atomJ->getType() != CollisionModel::Atom::AtomType::COM){

                        currentCharge = atomI->getCharge();

                    }else if (int(ceil(fabs(atomJ->getCharge()/Core::ELEMENTARY_CHARGE))) != 0 &&
                            atomI->getType() != CollisionModel::Atom::AtomType::COM){

                        currentCharge = atomJ->getCharge();

                    }
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
            }
            



            // Third contribution: ion <-> permanent dipole potential
            // This requires an ion and a dipole to be present
            // double dipoleDistanceScalar = 0;
            //double dipoleX = 0, dipoleY = 0, dipoleZ = 0;
            if(potentialsFF_ == "ALL" || potentialsFF_ == "PDIPOLE" || potentialsFF_ == "VDWPDIPOLE"){

                double currentCharge = 0;
                if(int(atomI->getCharge()/Core::ELEMENTARY_CHARGE) != 0 &&
                        moleculesPtr[1]->getIsDipole() == true){

                    currentCharge = atomI->getCharge();
                    /*dipoleX = moleculesPtr[1]->getDipole().x();
                    dipoleY = moleculesPtr[1]->getDipole().y();
                    dipoleZ = moleculesPtr[1]->getDipole().z();
                    dipoleDistanceScalar =  dipoleX * distance.x() +
                            dipoleY * distance.y() +
                            dipoleZ * distance.z();*/

                }else if (moleculesPtr[0]->getIsDipole() == true &&
                        int(atomJ->getCharge()/Core::ELEMENTARY_CHARGE) != 0){

                    currentCharge = atomJ->getCharge();
                    /*dipoleX = moleculesPtr[0]->getDipole().x();
                    dipoleY = moleculesPtr[0]->getDipole().y();
                    dipoleZ = moleculesPtr[0]->getDipole().z();
                    dipoleDistanceScalar =  dipoleX * distance.x() +
                            dipoleY * distance.y() +
                            dipoleZ * distance.z();*/
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
            }
            

            // Fourth contribution: quadrupole moment if background gas is N2
            // This requires an ion and N2 to be present
            if(potentialsFF_ == "ALL" || potentialsFF_ == "QUAD" || potentialsFF_ == "VDWQUAD" || potentialsFF_ == "VDWIIQUAD"){
                double partialChargeN2 = 0;
                double currentCharge = 0;
                if(int(ceil(fabs(atomI->getCharge()/Core::ELEMENTARY_CHARGE))) != 0 &&
                        (isN2 == true || isN2Approx == true)){

                    currentCharge = atomI->getCharge();
                    partialChargeN2 = atomJ->getPartCharge();

                }else if ((isN2 == true || isN2Approx == true) &&
                        int(ceil(fabs(atomJ->getCharge()/Core::ELEMENTARY_CHARGE))) != 0){

                    currentCharge = atomJ->getCharge();
                    partialChargeN2 = atomI->getPartCharge();
                }
                Core::Vector quadrupoleForce;
                quadrupoleForce.x(currentCharge * partialChargeN2 * 1./Core::ELECTRIC_CONSTANT * distance.x() / distanceCubed );
                quadrupoleForce.y(currentCharge * partialChargeN2 * 1./Core::ELECTRIC_CONSTANT * distance.y() / distanceCubed );
                quadrupoleForce.z(currentCharge * partialChargeN2 * 1./Core::ELECTRIC_CONSTANT * distance.z() / distanceCubed );
                forceMolecules[0] += quadrupoleForce;
                forceMolecules[1] += quadrupoleForce * (-1);
            }
        }
    }

    // add the C4 ion-induced dipole force contribution
    if(potentialsFF_ == "ALL" || potentialsFF_ == "II" || potentialsFF_ == "VDWII"){
        Core::Vector ionInducedForce;
        if(isN2Approx){
            ionInducedForce.x(1./(Core::ELECTRIC_CONSTANT) * collisionGasPolarizability_m3_/2 *
                    (eField[0]*eFieldDerivative[0] + eField[1]*eFieldDerivative[1] + eField[2]*eFieldDerivative[5]));
            ionInducedForce.y(1./(Core::ELECTRIC_CONSTANT) * collisionGasPolarizability_m3_/2 *
                    (eField[0]*eFieldDerivative[1] + eField[1]*eFieldDerivative[2] + eField[2]*eFieldDerivative[3]));
            ionInducedForce.z(1./(Core::ELECTRIC_CONSTANT) * collisionGasPolarizability_m3_/2 *
                    (eField[0]*eFieldDerivative[5] + eField[1]*eFieldDerivative[3] + eField[2]*eFieldDerivative[4]));
        }else{
            ionInducedForce.x(1./(Core::ELECTRIC_CONSTANT) * collisionGasPolarizability_m3_ *
                    (eField[0]*eFieldDerivative[0] + eField[1]*eFieldDerivative[1] + eField[2]*eFieldDerivative[5]));
            ionInducedForce.y(1./(Core::ELECTRIC_CONSTANT) * collisionGasPolarizability_m3_ *
                    (eField[0]*eFieldDerivative[1] + eField[1]*eFieldDerivative[2] + eField[2]*eFieldDerivative[3]));
            ionInducedForce.z(1./(Core::ELECTRIC_CONSTANT) * collisionGasPolarizability_m3_ *
                    (eField[0]*eFieldDerivative[5] + eField[1]*eFieldDerivative[3] + eField[2]*eFieldDerivative[4]));
        }
        forceMolecules[0] += ionInducedForce;
        forceMolecules[1] += ionInducedForce * (-1);
        forceII = ionInducedForce * (-1);
    }
}

double CollisionModel::MDForceField_Buckingham::calculateVDW(const CollisionModel::Atom& atomA, const CollisionModel::Atom& atomB, double distanceAbs){
    double distanceSquared = distanceAbs*distanceAbs;
    double distanceSquaredInverse = 1./distanceSquared;
    double distanceCubed = distanceSquared * sqrt(distanceSquared);
    
    double sigma = CollisionModel::Atom::calcLJSig(atomA, atomB);
    double sigma6 = sigma * sigma * sigma * sigma * sigma * sigma;
    double epsilon = CollisionModel::Atom::calcLJEps(atomA, atomB);
    double ljFactor = epsilon * (1.84e5 * exp(-12*distanceAbs/sigma) - 
                                2.25 * sigma6 * distanceSquaredInverse*distanceSquaredInverse*distanceSquaredInverse);

    return ljFactor;

}

void CollisionModel::MDForceField_Buckingham::populateInteractionTable(std::vector<Core::Particle*> particlesPtrs, 
                            std::unordered_map<std::string, std::shared_ptr<CollisionModel::MolecularStructure>> structureMap, 
                            std::string collisionGasIdentifier){
    for(const auto& particleI: particlesPtrs){
        for(const auto& atomI: particleI->getMolecularStructure()->getAtoms()){
            for(const auto& atomJ: (structureMap.at(collisionGasIdentifier))->getAtoms()){
                double sigma = CollisionModel::Atom::calcLJSig(*atomI, *atomJ);
                double ljFactor = Core::goldenSectionSearch<const CollisionModel::Atom&>(*atomI, *atomJ, calculateVDW, 1e-12, sigma, 1e-22);
                interactionMap.insert(*atomI, *atomJ, ljFactor);
                
            }
        }
    }
}

