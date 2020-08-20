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
 RS_VantHoffReaction.hpp

 A Van't Hoff type reaction:
 The reaction rate constant is estimated with the Van't Hoff reaction isobar from the reaction enthalpy,
 the equilibrium constant for the reaction at standard conditions and a reaction rate for the backward reaction

 ****************************/

#ifndef RS_VantHoffDependentReaction_hpp
#define RS_VantHoffDependentReaction_hpp

#include "RS_AbstractReaction.hpp"
#include "RS_constants.hpp"
#include <cmath>

namespace RS {
    class VantHoffReaction : public AbstractReaction {

    private:
        double H_R_;
        double K_s_;
        double k_backward_;

    public:
        VantHoffReaction(
                std::map<Substance*,int> educts,
                std::map<Substance*,int> products,
                double H_R,
                double K_s,
                double k_backward,
                std::string label
        );

        RS::ReactionEvent attemptReaction(ReactionConditions conditions, ReactiveParticle* particle, double dt) const;
        RS::ReactionEvent attemptReaction(CollisionConditions conditions, ReactiveParticle* particle) const;
    };
}

#endif //RS_VantHoffDependentReaction_hpp
