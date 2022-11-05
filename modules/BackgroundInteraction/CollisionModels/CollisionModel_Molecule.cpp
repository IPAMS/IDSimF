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
#include <iostream>


/**
 * Creates an empty molecule with a center-of-mass position and velocity
 * @param comPos the center-of-mass position of the created molecule
 * @param comVel the center-of-mass velocity of the created molecule
 */
CollisionModel::Molecule::Molecule(const Core::Vector &comPos, const Core::Vector &comVel):
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
 * @param diam molecule diameter in m
 */
CollisionModel::Molecule::Molecule(const Core::Vector &comPos, const Core::Vector &comVel, 
                                    const Core::Vector &agls, std::vector<std::shared_ptr<CollisionModel::Atom>> &atms,
                                    double diam):
    centerOfMassPos(comPos),
    centerOfMassVel(comVel),
    angles(agls),
    atoms(atms),
    diameter(diam)
{
    this->calcDipole();
    this->calcMass();
    this->setIsDipole();
    this->setIsIon(); 
    this->rotateMolecule();
    this->atomCount = atoms.size();  
}

/**
 * @brief Creates a molecule given by the molecular structure. 
 *          The mass and dipole is calculated through the information given by the atoms and does not need 
 *          to be specified.
 * 
 * @param comPos the center-of-mass position of the created molecule
 * @param comVel the center-of-mass velocity of the created molecule
 * @param structure molecular structure
 */

CollisionModel::Molecule::Molecule(const Core::Vector &comPos, const Core::Vector &comVel, 
                                    std::shared_ptr<CollisionModel::MolecularStructure> structure):
    centerOfMassPos(comPos),
    centerOfMassVel(comVel),
    isDipole(structure->getIsDipole()),
    isIon(structure->getIsIon()),
    mass(structure->getMass()),
    dipole(structure->getDipole()),
    dipoleMag(structure->getDipoleMag()),
    atomCount(structure->getAtoms().size()),
    diameter(structure->getDiameter()), 
    molecularStructureName(structure->getName())
{
    atoms.resize(atomCount);
    for(size_t i = 0; i < atomCount; i++) {
        this->atoms.at(i) = std::make_shared<Atom>(*(structure->getAtoms().at(i)));
    }

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
    this->rotateMolecule();
}

/**
 * Sets new diameter of the molecule
 */
void CollisionModel::Molecule::setDiameter(double diam){
    this->diameter = diam;
}

/**
 * Sets name of the molecular structure 
 */
void CollisionModel::Molecule::setMolecularStructureName(std::string name){
    this->molecularStructureName = name;
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
 * Gets the vector of atoms belonging to the molecule
 */
std::vector<std::shared_ptr<CollisionModel::Atom>>& CollisionModel::Molecule::getAtoms(){
    return atoms;
}

/**
 * Gets the diameter
 */
double CollisionModel::Molecule::getDiameter() const{
    return diameter;
}

/**
 * Gets the identifier of the molecular structure 
 */
std::string CollisionModel::Molecule::getMolecularStructureName() const{
    return molecularStructureName;
}

/**
 * Calculates the current mass of the molecule
 */
void CollisionModel::Molecule::calcMass(){
    this->mass = 0;
    for(auto& atom : atoms){
        //std::cout << atom->getMass() << " " << (int)atom->getType() << std::endl;
        this->mass += atom->getMass();
    }
    
}

/**
 * Calculates the current dipole vector
 */
void CollisionModel::Molecule::calcDipole(){
    
    this->dipole = Core::Vector(0.0, 0.0, 0.0);
    for(auto& atom : atoms){
        Core::Vector relChargePos = atom->getRelativePosition() * atom->getPartCharge();
        this->dipole += relChargePos;
    }
    this->dipoleMag = this->dipole.magnitude();
    
}

/**
 * Adds an additional atom to the molecule. Mass and dipole are recalculated.
 */
void CollisionModel::Molecule::addAtom(std::shared_ptr<CollisionModel::Atom> atm){
    this->atoms.push_back(atm);
    this->atomCount++;
    this->calcDipole();
    this->calcMass();
    this->setIsDipole();
    this->setIsIon();   
}

/**
 * Removes an atom from the molecule. Mass and dipole are recalculated.
 */
void CollisionModel::Molecule::removeAtom(std::shared_ptr<CollisionModel::Atom> atm){
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
 * Determines if the molecule is a dipole and sets the flag accordingly.
 */
void CollisionModel::Molecule::setIsDipole(){


    if(this->dipoleMag > 0 || this->dipoleMag < 0){
        this->isDipole = true;
    }else{
        this->isDipole = false;
    }
    
}

/**
 * Determines if the molecule is an ion and sets the flag accordingly.
 */
void CollisionModel::Molecule::setIsIon(){

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

/**
 * Rotates the molecule by an angle set beforehand 
 */
void CollisionModel::Molecule::rotateMolecule(){
    for(auto& atom : atoms){
       atom->rotate(this->angles);
    }
}