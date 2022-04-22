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

#include "CollisionModel_Molecule.hpp"
#include <algorithm>
#include <cmath>

/**
 * Creates an empty molecule with a center-of-mass position and velocity
 * @param comPos the center-of-mass position of the created molecule
 * @param comVel the center-of-mass velocity of the created molecule
 */
CollisionModel::Molecule::Molecule(Core::Vector &comPos, Core::Vector &comVel):
    centerOfMassPos(comPos),
    centerOfMassVel(comVel)
{}

/**
 * @brief Creates a molecule with given atoms and the starting angle configuration (rotation). 
 *          The mass and dipole is calculated through the information given by the atoms and does not need 
 *          to be specified.
 * 
 * @param comPos the center-of-mass position of the created molecule
 * @param comVel the center-of-mass velocity of the created molecule
 * @param agls the angles in x-y-z of the atoms relative position w.r.t the center-of-mass 
 * @param atms the vector with atoms contained in the molecule
 */
CollisionModel::Molecule::Molecule(Core::Vector &comPos, Core::Vector &comVel, 
                                    Core::Vector &agls, std::vector<CollisionModel::Atom*> atms):
    centerOfMassPos(comPos),
    centerOfMassVel(comVel),
    angles(agls),
    atoms(atms)
{
    this->calcDipole();
    this->calcMass();
    this->setIsDipole();
    this->setIsIon();   
}

/**
 * Sets new center-of-mass position
 */
void CollisionModel::Molecule::setComPos(Core::Vector comPos){
    this->centerOfMassPos = comPos;
}

/**
 * Sets new center-of-mass velocity
 */
void CollisionModel::Molecule::setComVel(Core::Vector comVel){
    this->centerOfMassVel = comVel;
}

/**
 * Sets new angles of the atoms w.r.t. the center-of-mass of the molecule
 */
void CollisionModel::Molecule::setAngles(Core::Vector agls){
    this->angles = agls;
}

/**
 * Gets the center-of-mass position
 */
Core::Vector& CollisionModel::Molecule::getComPos(){
    return centerOfMassPos;
}

/**
 * Gets the center-of-mass velocity
 */
Core::Vector& CollisionModel::Molecule::getComVel(){
    return centerOfMassVel;
}

/**
 * Gets the angles of the atoms w.r.t. the center-of-mass of the molecule
 */
Core::Vector& CollisionModel::Molecule::getAngles(){
    return angles;
}

/**
 * Returns if the molecule is a permanent dipole or not
 */
bool CollisionModel::Molecule::getIsDipole() const{
    return isDipole;
}

/**
 * Returns if the molecule is an ion or not
 */
bool CollisionModel::Molecule::getIsIon() const{
    return isIon;
}

/**
 * Gets mass of the molecule
 */
double CollisionModel::Molecule::getMass() const{
    return mass;
}

/**
 * Gets the dipole vector
 */
Core::Vector CollisionModel::Molecule::getDipole() const{
    return dipole;
}

/**
 * Gets the dipole magnitude
 */
double CollisionModel::Molecule::getDipoleMag() const{
    return dipoleMag;
}

/**
 * Gets the number of atoms belonging to the molecule
 */
std::size_t CollisionModel::Molecule::getAtomCount() const{
    return atomCount;
}

/**
 * Calculates the current mass of the molecule
 */
void CollisionModel::Molecule::calcMass(){
    
    for(auto* atom : atoms){
        this->mass += atom->getMass();
    }
    
}

/**
 * Calculates the current dipole vector
 */
void CollisionModel::Molecule::calcDipole(){
    
    for(auto* atom : atoms){
        Core::Vector relChargePos = atom->getRelativePosition() * atom->getPartCharge();
        this->dipole += relChargePos;
    }
    this->dipoleMag = this->dipole.magnitude();
    
}

/**
 * Adds an additional atom to the molecule. Mass and dipole are recalculated.
 */
void CollisionModel::Molecule::addAtom(CollisionModel::Atom* atm){
    this->atoms.push_back(atm);
    this->calcDipole();
    this->calcMass();
    this->setIsDipole();
    this->setIsIon();   
}

/**
 * Removes an atom from the molecule. Mass and dipole are recalculated.
 */
void CollisionModel::Molecule::removeAtom(CollisionModel::Atom* atm){
    auto atm_it = std::find(std::begin(atoms), std::end(atoms), atm);
    if(atm_it != std::end(atoms)) 
        atoms.erase(atm_it);
    this->calcDipole();
    this->calcMass();
    this->setIsDipole();
    this->setIsIon(); 
}

/**
 * Determines if the molecule is a dipole and sets the flag accordingly.
 */
void CollisionModel::Molecule::setIsDipole(){

    if(this->getIsIon() == 0){
        for(auto* atom : atoms){
            if(atom->getPartCharge() != 0){
                this->isDipole = 1;
                break;
            }
        }
        this->isDipole = 0;
    }else{
        this->isDipole = 0;
    }
    
}

/**
 * Determines if the molecule is an ion and sets the flag accordingly.
 */
void CollisionModel::Molecule::setIsIon(){

    double sumCharge = 0;
    for(auto* atom : atoms){
        sumCharge += atom->getCharge() / Core::ELEMENTARY_CHARGE;
    }
    double tol = 10E-15;
    if(sumCharge < tol && sumCharge > -tol){
        this->isIon = 0;
    }else{
        this->isIon = 1;
    }
    
}