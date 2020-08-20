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
 RS_FieldDependentVantHoffReaction.hpp

 Van't Hoff style reaction of an ion where the forward reaction rate is calculated from
 the effective ion temperature, estimated from ion mobility and electric field strength,
 and the reaction enthalpy, equilibrium constant between the forward and backward reaction at
 standard conditions and the reaction rate of the backward reaction.

 ****************************/

#ifndef RS_FieldDependentVantHoffReaction_hpp
#define RS_FieldDependentVantHoffReaction_hpp

#include "RS_AbstractReaction.hpp"
#include "RS_constants.hpp"
#include <cmath>

namespace RS {
    /**
     * @brief A field dependent van't Hoff reaction
     *
     * This class represents a van't Hoff style reaction of an ion where the forward reaction rate is calculated from
     * thermophysical parameters (reaction enthalpy and equilibrium constant) of the equilibrium between the
     * forward and backward reaction and the effective ion temperature.
     *
     * The effective ion temperature is calculated from the ion mobility (at standard conditions),
     * the molecular mass of the backgorund gas and the local background temperature and pressure.
     */
    class FieldDependentVantHoffReaction : public AbstractReaction {

    private:
        double H_R_;         ///< Reaction enthalphy for the forward reaction
        double K_s_;         ///< Equilibrium constant for the forward reaction
        double kBackward_;   ///< Rate constant for the backward reaction
        double mobility_;    ///< Electrical mobility of the charged particle
        double collisionGasMass_amu_; ///< Mass of the collision background gas in amu
        double collisionGasMass_kg_;  ///< Mass of the collision background gas in kg

        double const P0_pa_ = 101325; ///< The default pressure in Pa
        double const T0_K_ = 298.15;  ///< The default / standard temperature in K

    public:
        FieldDependentVantHoffReaction(
                std::map<Substance*,int> educts,
                std::map<Substance*,int> products,
                double H_R,
                double K_s,
                double kBackward,
                double electricMobility,
                double collisionGasMass,
                std::string label
        );

        RS::ReactionEvent attemptReaction(ReactionConditions conditions, ReactiveParticle* particle, double dt) const;
        RS::ReactionEvent attemptReaction(CollisionConditions conditions, ReactiveParticle* particle) const;
    };
}


#endif //RS_FieldDependentVantHoffReaction_hpp
