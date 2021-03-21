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

#include "RS_SimpleCollisionStepReaction.hpp"

RS::SimpleCollisionStepReaction::SimpleCollisionStepReaction(std::map<RS::Substance *, int> educts,
                                                             std::map<RS::Substance *, int> products,
                                                             double activationEnergy_eV, std::string label):
AbstractReaction(educts, products, true, "simple_step",label),
activationEnergy_(activationEnergy_eV / Core::JOULE_TO_EV)
{}

/*
 * This is a collision based reaction, thus the stochastic based probability is always zero
 * and this method should not be called
 */
RS::ReactionEvent RS::SimpleCollisionStepReaction::attemptReaction(RS::ReactionConditions conditions,
                                                                   ReactiveParticle *particle, double dt) const{
    throw ("Stochastic probability requested for collision based reaction SimpleCollisionStepReaction");
}

RS::ReactionEvent RS::SimpleCollisionStepReaction::attemptReaction(CollisionConditions conditions,
                                                                   ReactiveParticle *particle) const{
    if (conditions.totalCollisionEnergy > activationEnergy_) {
        return ReactionEvent{true, 1.0};
    } else {
        return ReactionEvent{false, 0.0};
    }
}

