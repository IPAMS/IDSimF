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
 RS_StaticThermalizingReaction.hpp

 Reaction which reinitializes the velocity of the simulated particle with a
 thermal velocity drawn from maxwell boltzmann distribution to model the effect of resonant charge transfer
 (Reactions of type A+ + A -> A + A+)

 The reaction probability is calculated from a reaction crossection

 ****************************/

#ifndef RS_CrossectionStaticReaction_hpp
#define RS_CrossectionStaticReaction_hpp

#include "RS_AbstractReaction.hpp"

namespace RS {
    class CrossectionThermalizingReaction : public AbstractReaction {

    public:
        CrossectionThermalizingReaction(
            const std::map<Substance*,int>& educts,
            const std::map<Substance*,int>& products,
            double reactionCrossectionM2,
            std::string label
        );

        //double rateConstant(ReactionConditions) const;
        RS::ReactionEvent attemptReaction(ReactionConditions conditions, ReactiveParticle* particle, double dt) const override;
        RS::ReactionEvent attemptReaction(CollisionConditions conditions, ReactiveParticle* particle) const override;

    private:
        const double PI_SQRT = std::sqrt(M_PI);
        double reactionCrossectionM2_ = 0.0;
        double reactionPartnerMassAmu_ = 0.0;
        double reactionPartnerMassKg_ = 0.0;
        double reactionPartnerNumberConcentration_1m3_ = 0.0;
    };
}


#endif //RS_ThermalizingStaticReaction_hpp
