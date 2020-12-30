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
 RS_AbstractReaction.hpp

 This implements an abstract chemical reaction in the RS simulation

 ****************************/


#ifndef RS_AbstractReaction_hpp
#define RS_AbstractReaction_hpp

#include "RS_Substance.hpp"
#include "RS_ReactiveParticle.hpp"
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

namespace RS{class AbstractReaction;}
std::ostream& operator<<(std::ostream& os, const RS::AbstractReaction& reac);

namespace RS {

    /**
     * Conditions for a chemical reaction (the parameters to calculate the actual
     * reaction rate and if a reaction actually occurs under defined conditions)
     */
    struct ReactionConditions {
        double temperature;   ///< the background temperature at the location of reaction (K)
        double electricField; ///< the background electric field at the location of reaction (V/m)
        double pressure;      ///< the background pressure at the location of reaction (Pa)
        double totalReactionEnergy; ///< total energy available in an reaction event (J)
    };

    /**
     * Conditions for an individual collision event
     */
    struct CollisionConditions {
        double totalCollisionEnergy; ///< total energy available in the collision (typically in the center of mass frame)
    };


    /**
     * Defining parameters for a reaction event.
     */
    struct ReactionEvent {
        bool reactionHappened;       ///< true if the reaction actually has taken place
        double reactionProbability;  ///< the total reaction probability of the individual reaction event
    };

    /**
     * An abstract chemical (elementary) reaction, base class for chemical reactions. It has educts and products and can calculate
     * somehow a reaction rate / rate constant (the details are implemented in the non virtual / non abstract
     * derived classes of this abstract base class).
     *
     * This class is intended to model elementary reactions, thus all stochiometric factors are integer values.
     *
     * The chemical substances are stored as pointers to Substances, thus all the comparisons and
     * equality checks in this class are based on the address identity of a substance.
     *
     */
    class AbstractReaction {
    private:
        std::map<Substance*,int> educts_;           ///< Map of all educts
        std::map<Substance*,int> products_;         ///< Map of all products
        std::map<Substance*,int> discreteEducts_;   ///< Map of all discrete educts (educts modeled as discrete particles)
        std::map<Substance*,int> discreteProducts_; ///< Map of all discrete products (products modeled as discrete particles)

        double staticReactionConcentration_; ///< The static reaction concentration (product of all static / isotropic educt concentrations)

        std::string typeLabel_; ///< A textual identifier for the type of the reaction
        std::string label_; ///< A textual label (which is also used for equality check of reactions)

        bool independent_;  ///< "Independent" reactions depend only on one discrete educt
        bool collisionReaction_; ///< Collision reactions are modeled by discrete collision events of a particle tracing model

    public:
        AbstractReaction(
                 std::map<Substance*,int> educts,
                 std::map<Substance*,int> products,
                 bool isCollisionReaction,
                 std::string typeLabel,
                 std::string label
                 );
        virtual ~AbstractReaction() = default;

        std::string getLabel() const;
        std::string getTypeLabel() const;
        //virtual double rateConstant(ReactionConditions) const = 0;
        bool generateRandomDecision(double probability) const;
        //a particle pointer is passed to the attemptReaction methods to allow modifications of the particles in the reaction
        virtual ReactionEvent attemptReaction(ReactionConditions conditions, ReactiveParticle* particle, double dt) const = 0;
        virtual ReactionEvent attemptReaction(CollisionConditions conditions, ReactiveParticle* particle) const = 0;

        bool isIndependent() const;
        bool isCollisionReaction() const;
        double staticReactionConcentration() const;
        void updateStaticReactionConcentration();
        std::map<Substance*,int> products() const;
        std::map<Substance*,int> educts() const;
        const std::map<Substance*,int>* discreteProducts() const;
        const std::map<Substance*,int>* discreteEducts() const;

        friend std::ostream& ::operator<<(std::ostream& os, const RS::AbstractReaction& reac);
    };
}

#endif /* RS_AbstractReaction_hpp */
