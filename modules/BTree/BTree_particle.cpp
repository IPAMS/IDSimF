/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2020 - Physical and Theoretical Chemistry /
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

#include "BTree_particle.hpp"
#include <iostream>

/**
 * Creates an empty particle
 */
BTree::Particle::Particle()
        :
        location_(Core::Vector(0.0, 0.0, 0.0)),
        velocity_(Core::Vector(0.0, 0.0, 0.0)),
        acceleration_(Core::Vector(0.0, 0.0, 0.0)),
        charge_(0.0),
        mobility_(3.5e-4),
        mass_(0.0),
        diameter_(0.0),
        STP_meanFreePath_(0.0),
        STP_meanThermalVelocity_(0.0),
        active_(true),
        invalid_(false),
        timeOfBirth_(0.0),
        splatTime_(0.0),
        hostNode_(nullptr){}

/**
 * Creates a massless charged particle at a location
 * @param location the location of the created particle
 * @param chargeElemCharges the charge of the particle (in units of elementary charges)
 */
BTree::Particle::Particle(const Core::Vector &location, double chargeElemCharges)
        :
        location_(location),
        velocity_(Core::Vector(0.0, 0.0, 0.0)),
        acceleration_(Core::Vector(0.0, 0.0, 0.0)),
        charge_(chargeElemCharges * Core::ELEMENTARY_CHARGE),
        mobility_(3.5e-4),
        mass_(0.0),
        diameter_(0.0),
        STP_meanFreePath_(0.0),
        STP_meanThermalVelocity_(0.0),
        active_(true),
        invalid_(false),
        timeOfBirth_(0.0),
        splatTime_(0.0),
        hostNode_(nullptr){}

/**
 * Creates a charged particle at a location
 *
 * @param location the location of the created particle
 * @param velocity the velocity of the created particle
 * @param chargeElemCharges the charge of the particle (in units of elementary charges)
 * @param massAMU the mass of the charged particle (in atomic mass units)
 */
BTree::Particle::Particle(const Core::Vector &location, const Core::Vector &velocity, double chargeElemCharges, double massAMU)
        :
        location_(location),
        velocity_(velocity),
        acceleration_(Core::Vector(0.0, 0.0, 0.0)),
        charge_(chargeElemCharges * Core::ELEMENTARY_CHARGE),
        mobility_(3.5e-4),
        mass_(massAMU * Core::AMU_TO_KG),
        diameter_(0.0),
        STP_meanFreePath_(0.0),
        STP_meanThermalVelocity_(0.0),
        active_(true),
        invalid_(false),
        timeOfBirth_(0.0),
        splatTime_(0.0),
        hostNode_(nullptr){}

/**
 * Creates a charged particle at a location
 * @param location the location of the created particle
 * @param velocity the velocity of the created particle
 * @param chargeElemCharges the charge of the particle (in units of elementary charges)
 * @param massAMU the mass of the charged particle (in atomic mass units)
 * @param timeOfBirth the time when this particle was created
 */
BTree::Particle::Particle(const Core::Vector &location, const Core::Vector &velocity, double chargeElemCharges, double massAMU,
                          double timeOfBirth)
        :
        location_(location),
        velocity_(velocity),
        acceleration_(Core::Vector(0.0, 0.0, 0.0)),
        charge_(chargeElemCharges * Core::ELEMENTARY_CHARGE),
        mobility_(3.5e-4),
        mass_(massAMU * Core::AMU_TO_KG),
        diameter_(0.0),
        STP_meanFreePath_(0.0),
        STP_meanThermalVelocity_(0.0),
        active_(true),
        invalid_(false),
        timeOfBirth_(timeOfBirth),
        splatTime_(0.0),
        hostNode_(nullptr){}

/**
* Creates a charged particle at a location
* @param location the location of the created particle
* @param velocity the velocity of the created particle
* @param chargeElemCharges the charge of the particle (in units of elementary charges)
* @param massAMU the mass of the charged particle (in atomic mass units)
* @param collisionDiameterM the effective collision diameter of the particle (in m)
* @param timeOfBirth the time when this particle was created
*/
BTree::Particle::Particle(const Core::Vector &location, const Core::Vector &velocity, double chargeElemCharges, double massAMU,
        double collisionDiameterM, double timeOfBirth)
        :
        location_(location),
        velocity_(velocity),
        acceleration_(Core::Vector(0.0, 0.0, 0.0)),
        charge_(chargeElemCharges * Core::ELEMENTARY_CHARGE),
        mobility_(3.5e-4),
        mass_(massAMU * Core::AMU_TO_KG),
        diameter_(collisionDiameterM),
        STP_meanFreePath_(0.0),
        STP_meanThermalVelocity_(0.0),
        active_(true),
        invalid_(false),
        timeOfBirth_(timeOfBirth),
        splatTime_(0.0),
        hostNode_(nullptr){}


/**
 * Sets a new location
 */
void BTree::Particle::setLocation(const Core::Vector &location){
        this->location_ = location;
}

/**
 * Gets the location
 */
Core::Vector& BTree::Particle::getLocation(){
    return(location_);
}

/**
 * Sets a velocity
 */
void BTree::Particle::setVelocity(const Core::Vector &velocity){
    this->velocity_ = velocity;
}

/**
 * Gets the velocity
 */
Core::Vector& BTree::Particle::getVelocity(){
    return(velocity_);
}

/**
 * Sets the acceleration
 */
void BTree::Particle::setAcceleration(const Core::Vector &acceleration) {
    this->acceleration_ = acceleration;
}

/**
 * Gets the acceleration
 */
Core::Vector& BTree::Particle::getAcceleration(){
    return(acceleration_);
}

/**
 * Sets a new host node for this particle
 * @param newHostNode link to a Core node which is the new host node
 */
void BTree::Particle::setHostNode(BTree::AbstractNode* newHostNode){
    this->hostNode_ = newHostNode;
}

/**
 * Gets the host node of this particle
 */
BTree::AbstractNode* BTree::Particle::getHostNode(){
    return (hostNode_);
}

/**
 * Sets the external index of the particle
 */
void BTree::Particle::setIndex(int index){
    index_ = index;
}

/**
 * Gets the external index of the particle
 */
int BTree::Particle::getIndex(){
    return (index_);
}


/**
 * Sets the charge in units of elementary charges
 * @param chargeElemCharges the new charge in units of elementary charges
 */
void BTree::Particle::setChargeElementary(double chargeElemCharges){
    this->charge_ = chargeElemCharges*Core::ELEMENTARY_CHARGE;
}

/**
 * Gets the charge of the particle
 * @return charge of the particle (in Coulomb)
 */
double BTree::Particle::getCharge() const{
    return(charge_);
}

/**
 * Sets active state of the particle
 */
void BTree::Particle::setActive(bool active){
    this->active_ = active;
}

/**
 * Gets the active state of the particle
 * @return true if the particle is active
 */
bool BTree::Particle::isActive() const{
    return(active_);
}

/**
 * Sets the invalid state of the particle
 */
void BTree::Particle::setInvalid(bool invalid){
    this->invalid_ = invalid;
}

/**
 * Gets the invalid state of the particle
 * @return true if the particle is invalid
 */
bool BTree::Particle::isInvalid() const{
    return(invalid_);
}

/**
 * Gets a floating point particle attribute
 * @param key a textual key which identifies the attribute
 * @return the value of the attribute with the given key
 */
double BTree::Particle::getFloatAttribute(const std::string& key) const{
    return attributesFloat_.at(key);
}

/**
 * Sets a floating point particle attribute
 * @param key a textual key which identifies the attribute
 * @param value the new value of the attribute
 */
void BTree::Particle::setFloatAttribute(const std::string& key, double value) {
    attributesFloat_[key] = value;
}

/**
 * Gets an integer particle attribute
 * @param key a textual key which identifies the attribute
 * @return the value of the attribute with the given key
 */
int BTree::Particle::getIntegerAttribute(const std::string& key) const{
    return attributesInteger_.at(key);
}

/**
 * Sets an integer particle attribute
 * @param key a textual key which identifies the attribute
 * @param value the new value of the attribute
 */
void BTree::Particle::setIntegerAttribute(const std::string& key, int value) {
    attributesInteger_[key] = value;
}

/**
 * Accesses the array of auxiliary parameters for collision models
 */
std::array<double,3>& BTree::Particle::getAuxCollisionParams() {
    return auxCollisionParams_;
}

/**
 * Sets the mobility
 */
void BTree::Particle::setMobility(double mobility){
    this->mobility_ = mobility;
}

/**
 * Gets the mobility
 */
double BTree::Particle::getMobility() const{
    return(mobility_);
}

/**
 * Sets the mean free path at standard temperature and standard pressure
 */
void BTree::Particle::setMeanFreePathSTP(double meanFreePathSTP) {
    this->STP_meanFreePath_ = meanFreePathSTP;
}

/**
 * Gets the mean free path at standard temperature and standard pressure
 */
double BTree::Particle::getMeanFreePathSTP() const{
    return(STP_meanFreePath_);
}

/**
 * Sets the mean thermal velocity at standard tempeature and standard pressure
 */
void BTree::Particle::setMeanThermalVelocitySTP(double meanVelocitySTP) {
    this->STP_meanThermalVelocity_ = meanVelocitySTP;
}

/**
 * Gets the mean velocity at standard temperature and standard pressure
 */
double BTree::Particle::getMeanThermalVelocitySTP() const{
    return(STP_meanThermalVelocity_);
}

/**
 * Sets the mass
 * @param massAMU mass in atomic mass units
 */
void BTree::Particle::setMassAMU(double massAMU){
    this->mass_ = massAMU*Core::AMU_TO_KG;
}

/**
 * Gets the mass (in kg)
 * @return the mass in kg
 */
double BTree::Particle::getMass() const{
    return(mass_);
}

/**
 * Sets the particle diameter
 */
void BTree::Particle::setDiameter(double diameter) {
    this->diameter_ = diameter;
}

/**
 * Gets the particle diameter
 */
double BTree::Particle::getDiameter() const{
    return(diameter_);
}

/**
 * Sets the time of particle creation / "birth"
 */
void BTree::Particle::setTimeOfBirth(double timeOfBirth){
    this->timeOfBirth_ = timeOfBirth;
}

/**
 * Gets the time of particle creation / "birth"
 */
double BTree::Particle::getTimeOfBirth() const{
    return(timeOfBirth_);
}

/**
 * Sets the "splat" / termination time
 */
void BTree::Particle::setSplatTime(double splatTime){
    this->splatTime_ = splatTime;
}

/**
 * Gets the "splat" / termination time
 */
double BTree::Particle::getSplatTime() const{
    return (splatTime_);
}