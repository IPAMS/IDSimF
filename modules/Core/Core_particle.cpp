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

#include "Core_particle.hpp"
#include <iostream>

// // FIXME: copy and move functionality 
// Core::Particle& Core::Particle::operator=( const Core::Particle& part)
// {
//     if(&part == this){
//         return *this;
//     }
//     index_ = part.index_;
//     location_ = part.location_;
//     velocity_ = part.velocity_;
//     acceleration_ = part.acceleration_;
//     charge_ = part.charge_;
//     mobility_ = part.mobility_;
//     mass_ = part.mass_;
//     diameter_ = part.diameter_;
//     STP_meanFreePath_ = part.STP_meanFreePath_;
//     STP_meanThermalVelocity_ = part.STP_meanThermalVelocity_;
//     active_ = part.active_;
//     invalid_ = part.invalid_;
//     timeOfBirth_ = part.timeOfBirth_;
//     splatTime_ = part.splatTime_;
//     attributesFloat_ = part.attributesFloat_;
//     attributesInteger_ = part.attributesInteger_;
//     auxCollisionParams_ = part.auxCollisionParams_;
//     molstrPtr.reset( new CollisionModel::MolecularStructure( *part.molstrPtr ) );
//     return *this;
// }

// Core::Particle& Core::Particle::operator=( Core::Particle&& part)
// {   
//     if(&part == this){
//         return *this;
//     }
//     index_ = part.index_;
//     part.index_ = 0;
//     location_ = part.location_;
//     part.location_ = Core::Vector(0.0, 0.0, 0.0);
//     velocity_ = part.velocity_;
//     part.velocity_ = Core::Vector(0.0, 0.0, 0.0);
//     acceleration_ = part.acceleration_;
//     part.acceleration_ = Core::Vector(0.0, 0.0, 0.0);
//     charge_ = part.charge_;
//     part.charge_ = 0;
//     mobility_ = part.mobility_;
//     part.mobility_ = 0;
//     mass_ = part.mass_;
//     part.mass_ = 0;
//     diameter_ = part.diameter_;
//     part.diameter_ = 0;
//     STP_meanFreePath_ = part.STP_meanFreePath_;
//     part.STP_meanFreePath_ = 0;
//     STP_meanThermalVelocity_ = part.STP_meanThermalVelocity_;
//     part.STP_meanThermalVelocity_ = 0;
//     active_ = part.active_;
//     part.active_ = true;
//     invalid_ = part.invalid_;
//     part.invalid_ = false;
//     timeOfBirth_ = part.timeOfBirth_;
//     part.timeOfBirth_ = 0;
//     splatTime_ = part.splatTime_;
//     part.splatTime_ = 0;
//     attributesFloat_ = part.attributesFloat_;
//     part.attributesFloat_ = {};
//     attributesInteger_ = part.attributesInteger_;
//     part.attributesInteger_ = {};
//     auxCollisionParams_ = part.auxCollisionParams_;
//     part.auxCollisionParams_ = {0.0, 0.0, 0.0};
//     molstrPtr = std::move( part.molstrPtr );
//     return *this;
// }

// Core::Particle::Particle(const Core::Particle& part) :
//     index_(part.index_),
//     location_(part.location_),
//     velocity_(part.velocity_),
//     acceleration_(part.acceleration_),
//     charge_(part.charge_),
//     mobility_(part.mobility_),
//     mass_(part.mass_),
//     diameter_(part.diameter_),
//     STP_meanFreePath_(part.STP_meanFreePath_),
//     STP_meanThermalVelocity_(part.STP_meanThermalVelocity_),
//     active_(part.active_),
//     invalid_(part.invalid_),
//     timeOfBirth_(part.timeOfBirth_),
//     splatTime_(part.splatTime_),
//     attributesFloat_(part.attributesFloat_),
//     attributesInteger_(part.attributesInteger_),
//     auxCollisionParams_(part.auxCollisionParams_),
//     molstrPtr(new CollisionModel::MolecularStructure(*part.molstrPtr))
// {}

// Core::Particle::Particle(Core::Particle&& part) :
//     index_(part.index_),
//     location_(part.location_),
//     velocity_(part.velocity_),
//     acceleration_(part.acceleration_),
//     charge_(part.charge_),
//     mobility_(part.mobility_),
//     mass_(part.mass_),
//     diameter_(part.diameter_),
//     STP_meanFreePath_(part.STP_meanFreePath_),
//     STP_meanThermalVelocity_(part.STP_meanThermalVelocity_),
//     active_(part.active_),
//     invalid_(part.invalid_),
//     timeOfBirth_(part.timeOfBirth_),
//     splatTime_(part.splatTime_),
//     attributesFloat_(part.attributesFloat_),
//     attributesInteger_(part.attributesInteger_),
//     auxCollisionParams_(part.auxCollisionParams_)
// {
//     part.index_ = 0;
//     part.location_ = Core::Vector(0.0, 0.0, 0.0);
//     part.velocity_ = Core::Vector(0.0, 0.0, 0.0);
//     part.acceleration_ = Core::Vector(0.0, 0.0, 0.0);
//     part.charge_ = 0;
//     part.mobility_ = 0;
//     part.mass_ = 0;
//     part.diameter_ = 0;
//     part.STP_meanFreePath_ = 0;
//     part.STP_meanThermalVelocity_ = 0;
//     part.active_ = true;
//     part.invalid_ = false;
//     part.timeOfBirth_ = 0;
//     part.splatTime_ = 0;
//     part.attributesFloat_ = {};
//     part.attributesInteger_ = {};
//     part.auxCollisionParams_ = {0.0, 0.0, 0.0};
//     molstrPtr = std::move( part.molstrPtr );
// }


/**
 * Creates a massless charged particle at a location
 * @param location the location of the created particle
 * @param chargeElemCharges the charge of the particle (in units of elementary charges)
 */
Core::Particle::Particle(const Core::Vector &location, double chargeElemCharges):
    location_(location),
    charge_(chargeElemCharges * Core::ELEMENTARY_CHARGE)
{}

/**
 * Creates a charged particle at a location
 *
 * @param location the location of the created particle
 * @param velocity the velocity of the created particle
 * @param chargeElemCharges the charge of the particle (in units of elementary charges)
 * @param massAMU the mass of the charged particle (in atomic mass units)
 */
Core::Particle::Particle(const Core::Vector &location, const Core::Vector &velocity,
                          double chargeElemCharges, double massAMU):
    location_(location),
    velocity_(velocity),
    charge_(chargeElemCharges * Core::ELEMENTARY_CHARGE),
    mass_(massAMU * Core::AMU_TO_KG)
{}

/**
 * Creates a charged particle at a location
 * @param location the location of the created particle
 * @param velocity the velocity of the created particle
 * @param chargeElemCharges the charge of the particle (in units of elementary charges)
 * @param massAMU the mass of the charged particle (in atomic mass units)
 * @param timeOfBirth the time when this particle was created
 */
Core::Particle::Particle(const Core::Vector &location, const Core::Vector &velocity, double chargeElemCharges,
                          double massAMU, double timeOfBirth):
    location_(location),
    velocity_(velocity),
    charge_(chargeElemCharges * Core::ELEMENTARY_CHARGE),
    mass_(massAMU * Core::AMU_TO_KG),
    timeOfBirth_(timeOfBirth)
{}

/**
* Creates a charged particle at a location
* @param location the location of the created particle
* @param velocity the velocity of the created particle
* @param chargeElemCharges the charge of the particle (in units of elementary charges)
* @param massAMU the mass of the charged particle (in atomic mass units)
* @param collisionDiameterM the effective collision diameter of the particle (in m)
* @param timeOfBirth the time when this particle was created
*/
Core::Particle::Particle(const Core::Vector &location, const Core::Vector &velocity,
                          double chargeElemCharges, double massAMU,
                          double collisionDiameterM, double timeOfBirth):
    location_(location),
    velocity_(velocity),
    charge_(chargeElemCharges * Core::ELEMENTARY_CHARGE),
    mobility_(3.5e-4),
    mass_(massAMU * Core::AMU_TO_KG),
    diameter_(collisionDiameterM),
    timeOfBirth_(timeOfBirth)
{}


/**
* Creates a charged particle at a location
* @param location the location of the created particle
* @param velocity the velocity of the created particle
* @param chargeElemCharges the charge of the particle (in units of elementary charges)
* @param massAMU the mass of the charged particle (in atomic mass units)
* @param collisionDiameterM the effective collision diameter of the particle (in m)
* @param timeOfBirth the time when this particle was created
* @param moleculeStructure pointer to the molecule which should be used in the MD Collision Model 
*/
Core::Particle::Particle(const Core::Vector &location, const Core::Vector &velocity,
                          double chargeElemCharges, double massAMU,
                          double collisionDiameterM, double timeOfBirth, 
                          std::shared_ptr<CollisionModel::MolecularStructure> moleculeStructure):
    location_(location),
    velocity_(velocity),
    charge_(chargeElemCharges * Core::ELEMENTARY_CHARGE),
    mobility_(3.5e-4),
    mass_(massAMU * Core::AMU_TO_KG),
    diameter_(collisionDiameterM),
    timeOfBirth_(timeOfBirth),
    molstrPtr(moleculeStructure)
{}


/**
 * Sets a new location
 */
void Core::Particle::setLocation(Core::Vector location){
        this->location_ = location;
}

/**
 * Gets the location
 */
Core::Vector& Core::Particle::getLocation(){
    return(location_);
}

/**
 * Sets a velocity
 */
void Core::Particle::setVelocity(Core::Vector velocity){
    this->velocity_ = velocity;
}

/**
 * Gets the velocity
 */
Core::Vector& Core::Particle::getVelocity(){
    return(velocity_);
}

/**
 * Sets the acceleration
 */
void Core::Particle::setAcceleration(Core::Vector acceleration) {
    this->acceleration_ = acceleration;
}

/**
 * Gets the acceleration
 */
Core::Vector& Core::Particle::getAcceleration(){
    return(acceleration_);
}

/**
 * Sets the external index of the particle
 */
void Core::Particle::setIndex(size_t index){
    index_ = index;
}

/**
 * Gets the external index of the particle
 */
std::size_t Core::Particle::getIndex() const{
    return (index_);
}


/**
 * Sets the charge in units of elementary charges
 * @param chargeElemCharges the new charge in units of elementary charges
 */
void Core::Particle::setChargeElementary(double chargeElemCharges){
    this->charge_ = chargeElemCharges*Core::ELEMENTARY_CHARGE;
}

/**
 * Gets the charge of the particle
 * @return charge of the particle (in Coulomb)
 */
double Core::Particle::getCharge() const{
    return(charge_);
}

/**
 * Sets active state of the particle
 */
void Core::Particle::setActive(bool active){
    this->active_ = active;
}

/**
 * Gets the active state of the particle
 * @return true if the particle is active
 */
bool Core::Particle::isActive() const{
    return(active_);
}

/**
 * Sets the invalid state of the particle
 */
void Core::Particle::setInvalid(bool invalid){
    this->invalid_ = invalid;
}

/**
 * Gets the invalid state of the particle
 * @return true if the particle is invalid
 */
bool Core::Particle::isInvalid() const{
    return(invalid_);
}

/**
 * Gets a floating point particle attribute
 * @param key a textual key which identifies the attribute
 * @return the value of the attribute with the given key
 */
double Core::Particle::getFloatAttribute(const std::string& key) const{
    return attributesFloat_.at(key);
}

/**
 * Sets a floating point particle attribute
 * @param key a textual key which identifies the attribute
 * @param value the new value of the attribute
 */
void Core::Particle::setFloatAttribute(const std::string& key, double value) {
    attributesFloat_[key] = value;
}

/**
 * Gets an integer particle attribute
 * @param key a textual key which identifies the attribute
 * @return the value of the attribute with the given key
 */
int Core::Particle::getIntegerAttribute(const std::string& key) const{
    return attributesInteger_.at(key);
}

/**
 * Sets an integer particle attribute
 * @param key a textual key which identifies the attribute
 * @param value the new value of the attribute
 */
void Core::Particle::setIntegerAttribute(const std::string& key, int value) {
    attributesInteger_[key] = value;
}

/**
 * Accesses the array of auxiliary parameters for collision models
 */
std::array<double,3>& Core::Particle::getAuxCollisionParams() {
    return auxCollisionParams_;
}

/**
 * Sets the mobility
 */
void Core::Particle::setMobility(double mobility){
    this->mobility_ = mobility;
}

/**
 * Gets the mobility
 */
double Core::Particle::getMobility() const{
    return(mobility_);
}

/**
 * Sets the mean free path at standard temperature and standard pressure
 */
void Core::Particle::setMeanFreePathSTP(double meanFreePathSTP) {
    this->STP_meanFreePath_ = meanFreePathSTP;
}

/**
 * Gets the mean free path at standard temperature and standard pressure
 */
double Core::Particle::getMeanFreePathSTP() const{
    return(STP_meanFreePath_);
}

/**
 * Sets the mean thermal velocity at standard tempeature and standard pressure
 */
void Core::Particle::setMeanThermalVelocitySTP(double meanVelocitySTP) {
    this->STP_meanThermalVelocity_ = meanVelocitySTP;
}

/**
 * Gets the mean velocity at standard temperature and standard pressure
 */
double Core::Particle::getMeanThermalVelocitySTP() const{
    return(STP_meanThermalVelocity_);
}

/**
 * Sets the mass
 * @param massAMU mass in atomic mass units
 */
void Core::Particle::setMassAMU(double massAMU){
    this->mass_ = massAMU*Core::AMU_TO_KG;
}

/**
 * Gets the mass (in kg)
 * @return the mass in kg
 */
double Core::Particle::getMass() const{
    return(mass_);
}

/**
 * Sets the particle diameter
 */
void Core::Particle::setDiameter(double diameter) {
    this->diameter_ = diameter;
}

/**
 * Gets the particle diameter
 */
double Core::Particle::getDiameter() const{
    return(diameter_);
}

/**
 * Sets the time of particle creation / "birth"
 */
void Core::Particle::setTimeOfBirth(double timeOfBirth){
    this->timeOfBirth_ = timeOfBirth;
}

/**
 * Gets the time of particle creation / "birth"
 */
double Core::Particle::getTimeOfBirth() const{
    return(timeOfBirth_);
}

/**
 * Sets the "splat" / termination time
 */
void Core::Particle::setSplatTime(double splatTime){
    this->splatTime_ = splatTime;
}

/**
 * Gets the "splat" / termination time
 */
double Core::Particle::getSplatTime() const{
    return (splatTime_);
}

/**
 * Sets the molecular structure pointer
 */
void Core::Particle::setMolecularStructure(std::shared_ptr<CollisionModel::MolecularStructure> molecularStructurePtr){
    this->molstrPtr = molecularStructurePtr;
}

/**
 * Gets the molecular structure pointer 
 */
std::shared_ptr<CollisionModel::MolecularStructure> Core::Particle::getMolecularStructure(){
    return molstrPtr;
}