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

#include "RS_StaticReaction.hpp"


RS::StaticReaction::StaticReaction(
        const std::map<RS::Substance*, int>& educts,
        const std::map<RS::Substance*, int>& products,
        double rateConstant,
        const std::string label):
AbstractReaction(educts, products,false, "static", label),
rateConstant_(rateConstant)
{}

RS::ReactionEvent RS::StaticReaction::attemptReaction(RS::ReactionConditions /*conditions*/,
                                                      RS::ReactiveParticle* /*particle*/, double dt) const{

    double reactionProbability = rateConstant_ * this->staticReactionConcentration() * dt;
    bool reactionHappened = generateRandomDecision(reactionProbability);

    return ReactionEvent{reactionHappened, reactionProbability};
}

/*
 * This is a purely stochastic reaction, thus the collision based probability is always zero
 * and this method should not be called
 */
RS::ReactionEvent RS::StaticReaction::attemptReaction(CollisionConditions /*conditions*/,
                                                      RS::ReactiveParticle* /*particle*/) const{
    throw std::logic_error(
        "Collision based reaction probability requested for purely stochastic reaction StaticReaction");
}
