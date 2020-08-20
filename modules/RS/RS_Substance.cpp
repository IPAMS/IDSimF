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

#include "RS_Substance.hpp"

RS::Substance::Substance(std::string name)
        :
        charge_(0.0),
        mass_(0.0),
        mobility_(0.0),
        collisionDiameter_(0.0),
        staticConcentration_(0.0),
        name_(name)
{}

RS::Substance::Substance(std::string name, RS::Substance::substanceType type)
:
RS::Substance(name)
{
    type_ = type;
}


RS::Substance::Substance(std::string name, std::string typeLabel) noexcept(false)
:
RS::Substance(name)
{
    RS::Substance::substanceType sType;
    if (typeLabel =="discrete"){
        type_= RS::Substance::substanceType::discrete;
    }
    else if (typeLabel =="isotropic"){
        type_= RS::Substance::substanceType::isotropic;
    }
    else if (typeLabel =="field"){
        type_= RS::Substance::substanceType::field;
    }
    else{
        throw RS::SubstanceException("illegal substance type");
    }
}


/**
 Get the name of the substance
 
 \returns the name of the substance
 */
std::string RS::Substance::name() const{
    return this->name_;
}

/**
 Get the type of the substance
 
 \returns the type of the substance (isotropic, discrete or field)
 */
const RS::Substance::substanceType RS::Substance::type() const{
    return type_;
}


/**
 Get the charge of the substance
 
 \returns the mass of the substance (in units of elementary charges)
 */
double RS::Substance::charge() const{
    return charge_;
}

/**
 Set the charge of the substance
 
 \param newCharge the charge (in units of elementary charges) to set for the substance
 */
void RS::Substance::charge(double newCharge){
    charge_ = newCharge;
}

/**
 Get the mass of the substance
 
 \returns the mass of the substance (in units of AMU)
 */
double RS::Substance::mass() const{
    return mass_;
}

/**
 Set the mass of the substance
 
 \param newMass the mass to set for the substance (in units of AMU)
 */
void RS::Substance::mass(double newMass){
    mass_ = newMass;
}

/**
 Get the static concentration of the substance
 
 \returns the static concentration of the substance
 */
double RS::Substance::staticConcentration() const{
    return staticConcentration_;
}

/**
 Set the static concentration of the substance
 
 \param newMass the static concentration to set for the substance
 */
void RS::Substance::staticConcentration(double newStaticConcentration){
    staticConcentration_ = newStaticConcentration;
}

/**
 Get the electrical mobility of the substance
 
 \returns the electrical mobility of the substance
 */
double RS::Substance::mobility() const{
    return mobility_;
}

/**
 Set the effective collision diameter of the substance
 
 \param newMobility the effective collision diameter to set
 */
void RS::Substance::collisionDiameter(double newCollisionDiameter){
    collisionDiameter_ = newCollisionDiameter;
}

/**
 Get the effective collision diameter of the substance

 \returns the effective collision diameter of the substance
 */
double RS::Substance::collisionDiameter() const{
    return collisionDiameter_;
}

/**
 Set the electrical mobility of the substance

 \param newMobility the electrical mobility to set
 */
void RS::Substance::mobility(double newMobility){
    mobility_ = newMobility;
}

bool RS::operator<(const RS::Substance a, const RS::Substance b){
    return a.name() < b.name();
}

/**
 * Output stream operator
 *
 * @param os
 * @param subst
 * @return An std::ostream (for chaining the stream operators)
 */
std::ostream& operator<<(std::ostream& os, const RS::Substance& subst)
{
    os <<"Substance Name: "<< subst.name() << " Substance Mass: "<<subst.mass();
    return os;
}

