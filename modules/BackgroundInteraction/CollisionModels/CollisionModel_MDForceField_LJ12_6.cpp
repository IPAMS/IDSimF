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

#include "CollisionModel_MDForceField_LJ12_6.hpp"
#include <array>

CollisionModel::MDForceField_LJ12_6::MDForceField_LJ12_6(double collisionGasPolarizability_m3):
    collisionGasPolarizability_m3_(collisionGasPolarizability_m3)
{}

void CollisionModel::MDForceField_LJ12_6::calculateForceField(std::vector<CollisionModel::Molecule*>& moleculesPtr,
                                                              std::vector<Core::Vector>& forceMolecules) {

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
            Core::Vector absPosAtomI = ion->getComPos() + atomI->getRelativePosition();
            Core::Vector absPosAtomJ = bgGas->getComPos() + atomJ->getRelativePosition();


            Core::Vector distance = absPosAtomI - absPosAtomJ;

            if(distance.magnitude() > 100e-10){
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

    // add the C4 ion-induced dipole force contribution
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
}