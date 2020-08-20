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
 RS_StaticReaction.hpp

 This implements a static chemical reaction in the RS simulation, a static reaction has a fixed reaction
 rate constant (k), independent of temperature / pressure / electric field

 ****************************/

#ifndef RS_StaticReaction_hpp
#define RS_StaticReaction_hpp

#include "RS_AbstractReaction.hpp"

namespace RS {
    class StaticReaction : public AbstractReaction {

    private:
        double rateConstant_;

    public:
        StaticReaction(
            std::map<Substance*,int> educts,
            std::map<Substance*,int> products,
            double rateConstant,
            std::string label
        );

        //double rateConstant(ReactionConditions) const;
        RS::ReactionEvent attemptReaction(ReactionConditions conditions, ReactiveParticle* particle, double dt) const;
        RS::ReactionEvent attemptReaction(CollisionConditions conditions, ReactiveParticle* particle) const;

    };
}


#endif //RS_StaticReaction_hpp
