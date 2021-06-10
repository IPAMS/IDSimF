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

 ------------
 RS_SimpleCollisionStepReaction.hpp

 Simple Activation Energy based reaction, primarily for collision based kinetcs models
 (If total reaction energy is below activation energy: No Reaction, if total reaction energy is above activation
 energy: reaction probability is 1)

 ****************************/
#ifndef RS_SimpleCollisionStepReaction_hpp
#define RS_SimpleCollisionStepReaction_hpp

#include "RS_AbstractReaction.hpp"

namespace RS {
    class SimpleCollisionStepReaction: public AbstractReaction {

    public:
        SimpleCollisionStepReaction(
                std::map<Substance*,int> educts,
                std::map<Substance*,int> products,
                double activationEnergy_eV,
                std::string label
        );

        RS::ReactionEvent attemptReaction(ReactionConditions conditions, ReactiveParticle* particle, double dt) const override;
        RS::ReactionEvent attemptReaction(CollisionConditions conditions, ReactiveParticle* particle) const override;

    private:
        double activationEnergy_= 0.0;
    };
}


#endif //RS_SimpleCollisionStepReaction_hpp
