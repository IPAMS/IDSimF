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
 * @param sig the Lennard-Jones parameter sigma of the atom (in units of meter)
 * @param eps the Lennard-Jones parameter epsilon of the atom (in units of Joule)
 */
CollisionModel::Atom::Atom(const Core::Vector &relPos, double massAMU, double chargeElemCharges, double partChargeElemCharges, 
                            CollisionModel::Atom::AtomType element, double sig, double eps):
    relativePosition(relPos),
    mass(massAMU * Core::AMU_TO_KG),
    charge(chargeElemCharges * Core::ELEMENTARY_CHARGE),
    partialCharge(partChargeElemCharges * Core::ELEMENTARY_CHARGE),
    type(element),
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
Core::Vector& CollisionModel::Atom::getRelativePosition(){
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
 * @brief Rotates atoms w.r.t its parent molecule based on given rotation angles in the
 *          x-y-z coordinate system. The implicitly given bond lengths are preserved.
 *          The angles are given as the accumalative angles in relation to the standard 
 *          configuration of the molecule (all three angles are zero), e.g.
 *          rotate([0,0,pi/2]) rotates by pi/2 about z from the zero configuration and 
 *          not by an additional pi/2 from the current configuration 
 * 
 * @param angles The current x-y-z rotation angles of the molecule
 */
void CollisionModel::Atom::rotate(const Core::Vector &angles){
    // NOTE: construct inverse rotation matrix by hand 
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
    return (atm1.getSigma() + atm2.getSigma()) / 2;
}