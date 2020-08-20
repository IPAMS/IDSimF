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

#include "RS_AbstractReaction.hpp"
#include "Core_randomGenerators.hpp"

/**
 * The default constructor of an abstract reaction
 *
 * @param educts a map of references / pointers to the educt Substances and their stociometric factors
 * @param products a map of references / pointers to the product substances and their stochiometric factors
 * @param label a textual label of this reaction (which is also used for equality checks)
 */
RS::AbstractReaction::AbstractReaction(
                       std::map<Substance*,int> educts,
                       std::map<Substance*,int> products,
                       bool isCollisionReaction,
                       std::string typeLabel,
                       std::string label)
:
educts_(educts),
products_(products),
collisionReaction_(isCollisionReaction),
typeLabel_(typeLabel),
label_(label)
{
    //init educts table
    //search discrete educts, calculate static reaction probability
    int nDiscrete = 0; //number of found discrete educt reaction partners
    for (auto it = educts_.begin(); it != educts_.end(); ++it) {
        RS::Substance* subst =it->first;
        const int sto_factor = it->second; //stochiometric factor
        if (subst->type() == Substance::substanceType::discrete){
            nDiscrete++;
            discreteEducts_[subst] = sto_factor;
        }
    }
    this->updateStaticReactionConcentration();

    //if we have only one discrete educt, the whole reaction is classified as independent
    if (nDiscrete == 1 ){
        independent_ = true;
    }
    
    //init products table
    for (auto it = products_.begin(); it != products_.end(); ++it) {
        //currently only discrete products are allowed by the config file parser
        //but the reaction implementation should be general
        RS::Substance* subst = it->first;
        const int sto_factor = it->second; //stochiometric factor
        if (subst->type() == Substance::substanceType::discrete){
            discreteProducts_[subst] = sto_factor;
        }
    }
}


/**
 * Gets the textual label from this reaction
 * @return the text label of this reaction
 */
std::string RS::AbstractReaction::getLabel() const {
    return label_;
}

/**
 * Gets the textual type label from this reaction which identifies the type of this reaction
 * @return the type text label of this reaction
 */
std::string RS::AbstractReaction::getTypeLabel() const {
    return typeLabel_;
}
/**
 * Generate a random decision: The result is true with the given probability
 *
 * @param probability the probabilty for the result to be true
 * @return the random decision
 */
bool RS::AbstractReaction::generateRandomDecision(double probability) const{
    double rndVal = Core::globalRandomGenerator->uniformRealRndValue();
    return (rndVal < probability);
}

/**
 * Independent reactions are dependent on only one discrete educt / substance modeled in terms of discrete particles
 * @return if this reaction is independent
 */
bool RS::AbstractReaction::isIndependent() const {
    return independent_;
}

/**
 * Collision reactions are reactions which are modeled by individual collision events of an external particle tracing
 * model which is not part of RS instead of an averaged approach where the reaction probability is calculated by the
 * timestep length.
 *
 * @return if this reactiion is a collision based reaction
 */
bool RS::AbstractReaction::isCollisionReaction() const{
    return collisionReaction_;
}

/**
 * The static probability is the product of all concentrations of isotropic / static educt species for this reaction.
 *
 * @return the static reaction probability of this reaction
 */
double RS::AbstractReaction::staticReactionConcentration() const {
    return staticReactionConcentration_;
}

/**
 * Updates / reacalculates the product of static reaction partner concentrations
 * of this reaction from the current concentrations of the isotropic / static educts.
 */
void RS::AbstractReaction::updateStaticReactionConcentration() {
    staticReactionConcentration_ = 1.0;
    for (auto it = educts_.begin(); it != educts_.end(); ++it) {
        RS::Substance* subst =it->first;
        const int sto_factor = it->second; //stochiometric factor
        if (subst->type() == Substance::substanceType::isotropic){
            staticReactionConcentration_ = staticReactionConcentration_ * ( std::pow(subst->staticConcentration(),sto_factor) );
        }

    }
}

/**
 * Gets the product substance map of this reaction.
 *
 * @return the map of products with their stociometric factors
 */
std::map<RS::Substance*,int> RS::AbstractReaction::products() const {
    return products_;
}

/**
 * Gets the educt substance map of this reaction.
 *
 * @return the map of educts with their stociometric factors
 */
std::map<RS::Substance*,int> RS::AbstractReaction::educts() const {
    return educts_;
}

/**
 * Gets the substance map of all discrete products.
 *
 * @return the map of discrete (modeled in terms of discrete particles) substances
 */
std::map<RS::Substance*,int> RS::AbstractReaction::discreteProducts() const {
    return discreteProducts_;
}

/**
 *
 * Gets the substance map of all discrete educts.
 *
 * @return the map of discrete (modeled in terms of discrete particles) educts
 */
std::map<RS::Substance*,int> RS::AbstractReaction::discreteEducts() const{
    return discreteEducts_;
}

/**
 * Output stream operator
 *
 * @param os an output stream
 * @param reac a chemical reaction
 * @return the modified output stream (for chaining)
 */
std::ostream& operator<<(std::ostream& os, const RS::AbstractReaction& reac)
{
    /*os <<"Reaction, educts: "<<std::endl;
    for (const auto &it : reac.educts_) {
        os << it.second << " " << it.first.name()<< std::endl;
    }

    os <<"Reaction, discrete Educts: "<<std::endl;
    for (const auto &it : reac.discreteEducts_) {
        os << it.second << " " << it.first.name()<< std::endl;
    }

    os <<"Reaction, products: "<<std::endl;
    for (const auto &it : reac.products_) {
        os << it.second << " " << it.first.name()<< std::endl;
    }

    os <<"Reaction, discrete Products: "<<std::endl;
    for (const auto &it : reac.discreteProducts_) {
        os << it.second << " " << it.first.name()<< std::endl;
    }*/

    os << reac.getLabel();
    return os;
}