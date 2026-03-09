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

#include "CollisionModel_Atom.hpp"
#include <cmath>

/**
 * Creates an atom with mass and charge at a position relative to the parent molecule
 * @param relPos the relative position of the created atom
 * @param massAMU the mass of the atom (in units of atomic mass)
 * @param chargeElemCharges the charge of the atom (in units of elementary charges)
 */
CollisionModel::Atom::Atom(const Core::Vector &relPos, double massAMU, double chargeElemCharges):
    relativePosition(relPos),
    mass(massAMU * Core::AMU_TO_KG),
    charge(chargeElemCharges * Core::ELEMENTARY_CHARGE)
{}

/**
 * Creates an atom with mass, charge and partial charges at a position relative to the parent molecule
 * @param relPos the relative position of the created atom
 * @param massAMU the mass of the atom (in units of atomic mass)
 * @param chargeElemCharges the charge of the atom (in units of elementary charges)
 * @param partChargeElemCharges the partial charge of the atom (in units of elementary charge)
 */
CollisionModel::Atom::Atom(const Core::Vector &relPos, double massAMU, double chargeElemCharges, double partChargeElemCharges):
    relativePosition(relPos),
    mass(massAMU * Core::AMU_TO_KG),
    charge(chargeElemCharges * Core::ELEMENTARY_CHARGE),
    partialCharge(partChargeElemCharges * Core::ELEMENTARY_CHARGE)
{}

/**
 * Creates an atom with mass, charge and partial charges at a position relative to the parent molecule
 * @param relPos the relative position of the created atom
 * @param massAMU the mass of the atom (in units of atomic mass)
 * @param chargeElemCharges the charge of the atom (in units of elementary charges)
 * @param partChargeElemCharges the partial charge of the atom (in units of elementary charge)
 * @param element the atom type
 * @param speciesIndex a unique numerical ID for the atomic species of this atom
 * @param sig the Lennard-Jones parameter sigma of the atom (in units of meter)
 * @param eps the Lennard-Jones parameter epsilon of the atom (in units of Joule)
 */
CollisionModel::Atom::Atom(const Core::Vector &relPos, double massAMU, double chargeElemCharges, double partChargeElemCharges, 
                            CollisionModel::Atom::AtomType element, std::size_t speciesIndex, double sig, double eps):
    relativePosition(relPos),
    mass(massAMU * Core::AMU_TO_KG),
    charge(chargeElemCharges * Core::ELEMENTARY_CHARGE),
    partialCharge(partChargeElemCharges * Core::ELEMENTARY_CHARGE),
    type(element),
    atomSpeciesIndex_(speciesIndex),
    sigma(sig),
    epsilon(eps)
{}

/**
 * Sets new relative position w.r.t the parent molecule 
 */
void CollisionModel::Atom::setRelativePosition(Core::Vector relPos){
    this->relativePosition = relPos;
}

/**
 * Sets new mass 
 * @param massAMU the new mass in atomic units 
 */
void CollisionModel::Atom::setMass(double massAMU){
    this->mass = massAMU * Core::AMU_TO_KG;
}

/**
 * Sets new atom type
 */
void CollisionModel::Atom::setType(CollisionModel::Atom::AtomType element){
    this->type = element;
}

/**
 * Sets the numeric species index
 */
void CollisionModel::Atom::setSpeciesIndex(std::size_t speciesIndex) {
    this->atomSpeciesIndex_ = speciesIndex;
}

/**
 * Sets new LJ parameter sigma
 */
void CollisionModel::Atom::setSigma(double sig){
    this->sigma = sig;
}

/**
 * Sets new LJ paramter epsilon
 */
void CollisionModel::Atom::setEpsilon(double eps){
    this->epsilon = eps;
}

/**
 * Sets new charge
 * @param chargeElemCharges the new charge in units of elementary charge 
 */
void CollisionModel::Atom::setCharge(double chargeElemCharges){
    this->charge = chargeElemCharges * Core::ELEMENTARY_CHARGE;
}

/**
 * Sets new partial charge
 * @param partChargeElemCharges the new partial charge in units of elementary charge 
 */
void CollisionModel::Atom::setPartCharge(double partChargeElemCharges){
    this->partialCharge = partChargeElemCharges * Core::ELEMENTARY_CHARGE;
}

/**
 * Get realtive position w.r.t its parent molecule
 */
Core::Vector& CollisionModel::Atom::getRelativePosition_(){
    return relativePosition;
}  

/**
 * Get atom mass
 */
double CollisionModel::Atom::getMass() const{
    return mass;
}

/**
 * Get atom type 
 */
CollisionModel::Atom::AtomType CollisionModel::Atom::getType() const{
    return type;
}

/**
 * Get numeric atomic species index
 */
std::size_t CollisionModel::Atom::getSpeciesIndex() const {
    return atomSpeciesIndex_;
}

/**
 * Get LJ parameter sigma 
 */
double CollisionModel::Atom::getSigma() const{
    return sigma;
}

/**
 * Get LJ parameter epsilon 
 */
double CollisionModel::Atom::getEpsilon() const{
    return epsilon;
}

/**
 * Get atom charge 
 */
double CollisionModel::Atom::getCharge() const{
    return charge;
}

/**
 * Get partial atom charge 
 */
double CollisionModel::Atom::getPartCharge() const{
    return partialCharge;
}

/**
 * Calculates the approximate LJ interaction parameter epsilon between two atoms 
 */
double CollisionModel::Atom::calcLJEps(const CollisionModel::Atom &atm1, const CollisionModel::Atom &atm2){
    return std::sqrt(atm1.getEpsilon() * atm2.getEpsilon());
}

/**
 * Calculates the approximate LJ interaction parameter sigma between two atoms 
 */
double CollisionModel::Atom::calcLJSig(const CollisionModel::Atom &atm1, const CollisionModel::Atom &atm2){
    // if(atm1.getType() == CollisionModel::Atom::AtomType::COM || atm2.getType() == CollisionModel::Atom::AtomType::COM){
    //     return 0.0; 
    // }else{
    return (atm1.getSigma() + atm2.getSigma()) / 2;
    // }
}

CollisionModel::Atom::AtomType CollisionModel::Atom::from_string(std::string str) {
    if(str == "C") {
        return CollisionModel::Atom::AtomType::C;
    } else if(str == "O") {
        return CollisionModel::Atom::AtomType::O;
    } else if(str == "N") {
        return CollisionModel::Atom::AtomType::N;
    } else if(str == "H") {
        return CollisionModel::Atom::AtomType::H;
    } else if(str == "He") {
        return CollisionModel::Atom::AtomType::He;
    } else if(str == "Ar") {
        return CollisionModel::Atom::AtomType::Ar;
    } else if(str == "Cl") {
        return CollisionModel::Atom::AtomType::Cl;
    } else if(str == "Li") {
        return CollisionModel::Atom::AtomType::Li;
    } else if(str == "Sn") {
        return CollisionModel::Atom::AtomType::Sn;
    } else if(str == "Na") {
        return CollisionModel::Atom::AtomType::Na;
    } else if(str == "COM") {
        return CollisionModel::Atom::AtomType::COM;
    } else {
        return CollisionModel::Atom::AtomType::UN;
        //throw std::invalid_argument("No such AtomType can be found.");
    }
}