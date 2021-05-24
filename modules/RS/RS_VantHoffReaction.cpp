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

#include "RS_VantHoffReaction.hpp"

RS::VantHoffReaction::VantHoffReaction(
        std::map<Substance*, int> educts,
        std::map<Substance*, int> products,
        double H_R,
        double K_s,
        double k_backward,
        std::string label):
AbstractReaction(educts, products, false, "vanthoff", label),
H_R_(H_R),
K_s_(K_s),
k_backward_(k_backward)
{}


RS::ReactionEvent RS::VantHoffReaction::attemptReaction(RS::ReactionConditions conditions,
                                                        ReactiveParticle* /*particle*/,
                                                        double dt) const{
    double k_forward =
            1.0 /
            (std::exp(H_R_ / RGas * (1.0 / conditions.temperature - 1.0 / T_standard)) * K_s_)
            * k_backward_;

    double reactionProbability = k_forward* this->staticReactionConcentration() * dt;
    bool reactionHappened = generateRandomDecision(reactionProbability);

    return ReactionEvent{reactionHappened, reactionProbability};
}

/*
 * This is a purely stochastic reaction, thus the collision based probability is always zero
 * and this method should not be called
 */
RS::ReactionEvent RS::VantHoffReaction::attemptReaction(
        CollisionConditions /*conditions*/, ReactiveParticle* /*particle*/) const{
    throw ("Collision based reaction probability requested for purely stochastic reaction VantHoffReaction");
}
