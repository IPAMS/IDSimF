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

#include "RS_CrossectionThermalizingReaction.hpp"
#include "RS_util.hpp"

RS::CrossectionThermalizingReaction::CrossectionThermalizingReaction(
        const std::map<RS::Substance*, int>& educts,
        const std::map<RS::Substance*, int>& products,
        double reactionDiamM,
        const std::string label):
AbstractReaction(educts, products, false, "thermalizing", label){
    reactionCrossectionM2_ = M_PI * (reactionDiamM/2.0) * (reactionDiamM/2.0);
    reactionPartnerNumberConcentration_1m3_ = staticReactionConcentration();

    // check if reaction has the correct structure of educts:
    //auto educts = this->educts();
    int nIsotropic = 0;
    reactionPartnerMassAmu_ = 0.0;

    for (auto it = educts.begin(); it != educts.end(); ++it) {
        RS::Substance* subst =it->first;
        if (subst->type() == Substance::substanceType::isotropic){
            nIsotropic++;
            reactionPartnerMassAmu_ = subst->mass();
        }
    }
    if (nIsotropic != 1) {
        throw std::invalid_argument("Cross section thermalizing reaction requires exactly one isotropic substance");
    }

    if (reactionPartnerMassAmu_ <= 0.0) {
        throw std::invalid_argument("Isotropic reaction partner requires molecular mass for cross section thermalizing reaction");
    }

    reactionPartnerMassKg_ = reactionPartnerMassAmu_ * Core::AMU_TO_KG;
}

RS::ReactionEvent RS::CrossectionThermalizingReaction::attemptReaction(RS::ReactionConditions conditions,
                                                                  RS::ReactiveParticle *particle, double dt) const{

    // calculate reaction probabilty from reaction partner density and reaction cross section
    // in a hard sphere way (analogous to collision probability in HS collision model):

    double vParticle = particle->getVelocity().magnitude(); //reactive particle velocity

    // Calculate the mean free path (MFP) from current particle velocity:

    // a static particle leads in static gas leads to a relative velocity of zero, which leads
    // to undefined behavior due to division by zero later.
    // The whole process converges to the MFP and collision probability of a static particle, thus
    // it is possible to assume a small velocity (1 nm/s) for the static ions to get rid of undefined behavior
    if (vParticle < 1e-9){
        vParticle = 1e-9;
    }

    // Calculate the mean reaction partner gas speed (m/s)
    double temperature_K = conditions.temperature;
    double vMeanGas = std::sqrt(8.0*Core::K_BOLTZMANN*temperature_K/M_PI/reactionPartnerMassKg_);

    // Calculate the median gas speed (m/s)
    double vMedianGas = std::sqrt(2.0*Core::K_BOLTZMANN*temperature_K/reactionPartnerMassKg_);

    // Compute the mean relative speed (m/s) between ion and gas.
    double s = vParticle / vMedianGas;
    double cMeanRel = vMeanGas * (
            (s + 1.0/(2.0*s)) * 0.5 * PI_SQRT * std::erf(s) + 0.5 * std::exp(-s*s) );

    // Compute mean-free-path (m)
    double effectiveMFP_m = (vParticle / cMeanRel) / (reactionPartnerNumberConcentration_1m3_ * reactionCrossectionM2_);

    // k*T*v/c / (p * sigma) == v/c * (k*T)/p * 1/sigma = v * 1/N *1/sigma = v/c/(N*sigma) = lambda
    // TODO: DImensional analysis

    // Compute probability of collision in the current time-step.
    double reactionProbability = 1.0 - std::exp(-vParticle* dt / effectiveMFP_m);
    bool reactionHappened = generateRandomDecision(reactionProbability);

    if (reactionHappened){
        //reinitalize the reacting particle with a random MB velocity

        //first: get the product particle mass (since the particle will be updated with the new chemical
        //species afterwards and outside of this method
        double productMass = this->discreteProducts()->begin()->first->mass();

        particle->setVelocity(
                RS::util::maxwellBoltzmannRandomVelocity(conditions.temperature,productMass));
    }

    return ReactionEvent{reactionHappened, reactionProbability};
}

/*
 * This is a purely stochastic reaction, thus the collision based probability is always zero
 * and this method should not be called
 */
RS::ReactionEvent RS::CrossectionThermalizingReaction::attemptReaction(CollisionConditions, RS::ReactiveParticle*) const{
    throw std::logic_error(
            "Collision based reaction probability requested for purely stochastic reaction Thermalizing Reaction");
}
