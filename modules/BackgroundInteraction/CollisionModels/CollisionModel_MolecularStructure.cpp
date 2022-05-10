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

#include "CollisionModel_MolecularStructure.hpp"
#include <algorithm>

/**
 * @brief Creates a molecular structur with given atoms. 
 *          The mass and dipole is calculated through the information given by the atoms and does not need 
 *          to be specified.
 * 
 * @param atms the vector with atoms contained in the molecule
 * @param diam molecule diameter in m
 */
CollisionModel::MolecularStructure::MolecularStructure(std::vector<std::shared_ptr<CollisionModel::Atom>> atms, double diam):
    atoms(atms),
    diameter(diam)
{
    this->calcDipole();
    this->calcMass();
    this->setIsDipole();
    this->setIsIon(); 
    this->atomCount = atoms.size();  
}

/**
 * Sets new diameter of the MolecularStructure
 */
void CollisionModel::MolecularStructure::setDiameter(double diam){
    this->diameter = diam;
}

/**
 * Returns if the MolecularStructure is a permanent dipole or not
 */
bool CollisionModel::MolecularStructure::getIsDipole() const{
    return isDipole;
}

/**
 * Returns if the MolecularStructure is an ion or not
 */
bool CollisionModel::MolecularStructure::getIsIon() const{
    return isIon;
}

/**
 * Gets mass of the MolecularStructure
 */
double CollisionModel::MolecularStructure::getMass() const{
    return mass;
}

/**
 * Gets the dipole vector
 */
Core::Vector CollisionModel::MolecularStructure::getDipole() const{
    return dipole;
}

/**
 * Gets the dipole magnitude
 */
double CollisionModel::MolecularStructure::getDipoleMag() const{
    return dipoleMag;
}

/**
 * Gets the number of atoms belonging to the MolecularStructure
 */
std::size_t CollisionModel::MolecularStructure::getAtomCount() const{
    return atomCount;
}

/**
 * Gets the vector of atoms belonging to the MolecularStructure
 */
std::vector<std::shared_ptr<CollisionModel::Atom>> CollisionModel::MolecularStructure::getAtoms() const{
    return atoms;
}

/**
 * Gets the diameter
 */
double CollisionModel::MolecularStructure::getDiameter() const{
    return diameter;
}

/**
 * Calculates the current mass of the MolecularStructure
 */
void CollisionModel::MolecularStructure::calcMass(){
    this->mass = 0;
    for(auto& atom : atoms){
        this->mass += atom->getMass();
    }
    
}

/**
 * Calculates the current dipole vector
 */
void CollisionModel::MolecularStructure::calcDipole(){
    
    this->dipole = Core::Vector(0.0, 0.0, 0.0);
    for(auto& atom : atoms){
        Core::Vector relChargePos = atom->getRelativePosition() * atom->getPartCharge();
        this->dipole += relChargePos;
    }
    this->dipoleMag = this->dipole.magnitude();
    
}

/**
 * Adds an additional atom to the MolecularStructure. Mass and dipole are recalculated.
 */
void CollisionModel::MolecularStructure::addAtom(std::shared_ptr<CollisionModel::Atom> atm){
    this->atoms.push_back(atm);
    this->atomCount++;
    this->calcDipole();
    this->calcMass();
    this->setIsDipole();
    this->setIsIon();   
}

/**
 * Removes an atom from the MolecularStructure. Mass and dipole are recalculated.
 */
void CollisionModel::MolecularStructure::removeAtom(std::shared_ptr<CollisionModel::Atom> atm){
    auto atm_it = std::find(std::begin(atoms), std::end(atoms), atm);
    if(atm_it != std::end(atoms)) 
        atoms.erase(atm_it);
    this->atomCount--;
    this->calcDipole();
    this->calcMass();
    this->setIsDipole();
    this->setIsIon(); 
}

/**
 * Determines if the MolecularStructure is a dipole and sets the flag accordingly.
 */
void CollisionModel::MolecularStructure::setIsDipole(){


    if(this->dipoleMag > 0 || this->dipoleMag < 0){
        this->isDipole = true;
    }else{
        this->isDipole = false;
    }
    
}

/**
 * Determines if the MolecularStructure is an ion and sets the flag accordingly.
 */
void CollisionModel::MolecularStructure::setIsIon(){

    double sumCharge = 0;
    for(auto& atom : atoms){
        sumCharge += atom->getCharge() / Core::ELEMENTARY_CHARGE;
    }
    double tol = 10E-15;
    if(sumCharge < tol && sumCharge > -tol){
        this->isIon = false;
    }else{
        this->isIon = true;
    }
    
}

// create an empty static collection of molecular structures
std::unordered_map<std::string, std::shared_ptr<CollisionModel::MolecularStructure> > CollisionModel::MolecularStructure::createCollection(){
    
    std::unordered_map<std::string, std::shared_ptr<CollisionModel::MolecularStructure> > molecularCollection;
    return molecularCollection;
}

std::unordered_map<std::string, std::shared_ptr<CollisionModel::MolecularStructure> > CollisionModel::MolecularStructure::molecularStructureCollection = CollisionModel::MolecularStructure::createCollection();