/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2022 - Physical and Theoretical Chemistry /
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

#include "CollisionModel_SoftSphere.hpp"
#include "Core_math.hpp"
#include "Core_randomGenerators.hpp"
#include "Core_utils.hpp"

/**
 * Constructor for static background gas pressure and temperature
 *
 * if maxwellianApproximation is true, the colliding background gas particle velocity is drawn
 * from a Maxwell Boltzmann distribution instead of the more correct relative velocity distribution
 * between charged particle and background gas
 */
CollisionModel::SoftSphereModel::SoftSphereModel(
        double staticPressure,
        double staticTemperature,
        double collisionGasMassAmu,
        double collisionGasDiameterM,
        bool maxwellianApproximation)
        :
        SoftSphereModel(
                getConstantDoubleFunction(staticPressure),
                getConstantVectorFunction(Core::Vector(0.0, 0.0, 0.0)),
                staticTemperature,
                collisionGasMassAmu,
                collisionGasDiameterM,
                maxwellianApproximation) { }

/**
* Constructor for static background gas pressure and temperature with a custom
* after collision function (most probably to model collision based reactions).
* After collision functions perform custom operations with the colliding particle and the
* conditions of the individual collision.
*
*/
CollisionModel::SoftSphereModel::SoftSphereModel(
        double staticPressure,
        double staticTemperature,
        double collisionGasMassAmu,
        double collisionGasDiameterM,
        std::function<void(RS::CollisionConditions, Core::Particle&)> afterCollisionFunction,
        bool maxwellianApproximation)
        :
        SoftSphereModel(
                getConstantDoubleFunction(staticPressure),
                getConstantVectorFunction(Core::Vector(0.0, 0.0, 0.0)),
                getConstantDoubleFunction(staticTemperature),
                collisionGasMassAmu,
                collisionGasDiameterM,
                std::move(afterCollisionFunction),
                maxwellianApproximation) { }

/**
 * Constructor for location dependent pressure and velocity, given as spatially resolved
 * functions.
 */
CollisionModel::SoftSphereModel::SoftSphereModel(
        std::function<double(Core::Vector& location)> pressureFunction,
        std::function<Core::Vector(Core::Vector& location)> velocityFunction,
        double staticTemperature,
        double collisionGasMassAmu,
        double collisionGasDiameterM,
        bool maxwellianApproximation)
        :
        SoftSphereModel(
                std::move(pressureFunction),
                std::move(velocityFunction),
                getConstantDoubleFunction(staticTemperature),
                collisionGasMassAmu,
                collisionGasDiameterM,
                nullptr,
                maxwellianApproximation) { }

/**
* Constructor for location dependent pressure and velocity, given as spatially resolved
* functions and custom after collision function (most probably to model collision based reactions).
*/
CollisionModel::SoftSphereModel::SoftSphereModel(
        std::function<double(Core::Vector& location)> pressureFunction,
        std::function<Core::Vector(Core::Vector& location)> velocityFunction,
        std::function<double(const Core::Vector&)>temperatureFunction,
        double collisionGasMassAmu,
        double collisionGasDiameterM,
        std::function<void(RS::CollisionConditions, Core::Particle&)> afterCollisionFunction,
        bool maxwellianApproximation)
        :
        maxwellianApproximation_(maxwellianApproximation),
        collisionGasMass_Amu_(collisionGasMassAmu),
        collisionGasMass_kg_(collisionGasMassAmu*Core::AMU_TO_KG),
        collisionGasDiameter_m_(collisionGasDiameterM),
        pressureFunction_(std::move(pressureFunction)),
        velocityFunction_(std::move(velocityFunction)),
        temperatureFunction_(std::move(temperatureFunction)),
        afterCollisionActionFunction_(std::move(afterCollisionFunction)) { }

void CollisionModel::SoftSphereModel::updateModelParticleParameters(Core::Particle& /*ion*/) const {}

void CollisionModel::SoftSphereModel::initializeModelParticleParameters(Core::Particle& /*ion*/) const {}

void CollisionModel::SoftSphereModel::updateModelTimestepParameters(int /*timestep*/, double /*time*/) {}

void CollisionModel::SoftSphereModel::modifyAcceleration(Core::Vector& /*acceleration*/, Core::Particle& /*ion*/,
                                                         double /*dt*/) {}

/**
 * Modify the velocity of a charged particle "ion" with a potential random collision in time step "dt"
 * (only one or none collision can happen, not multiple. Thus, the probability of multiple colllision events
 * happening in "dt" has to be low)
 */
void CollisionModel::SoftSphereModel::modifyVelocity(Core::Particle &ion, double dt) {

    Core::RandomSource *rndSource = Core::globalRandomGeneratorPool->getThreadRandomSource();

    Core::Vector pLocation = ion.getLocation();
    double localPressure_Pa = pressureFunction_(pLocation);

    if (Core::isDoubleEqual(localPressure_Pa, 0.0)) {
        return; //pressure 0 means no collision at all
    }

    // Transform the frame of reference in a frame where the mean background gas velocity is zero.
    Core::Vector vGasMean = velocityFunction_(pLocation);
    Core::Vector vFrameMeanBackRest = ion.getVelocity() - vGasMean;

    double vRelIonMeanBackRest = vFrameMeanBackRest.magnitude(); //relative ion relative to bulk gas velocity

    // Calculate collision cross section between particle and collision gas:
    // TODO: It seems to be unnecessary to constantly recalculate this
    // value, cache the calculated values somehow?
    // Correct the Sum of the Diameters variable according to kinetic energy
    double m1 = ion.getMass();
    double m2 = collisionGasMass_kg_;
    double reducedMass = (m1 * m2) / (m1 + m2);
    double kineticEnergy = 0.5 * reducedMass * (vRelIonMeanBackRest*vRelIonMeanBackRest);
    double vss_collision_omega = ion.getFloatAttribute(VSS_OMEGA);
    double DiameterCorrectionParameter = std::sqrt(std::pow(kineticEnergy, -vss_collision_omega));
    double sigma_m2;
    if(vss_collision_omega!=0.0){
        double correctedMeanDiameter = ((ion.getDiameter() + collisionGasDiameter_m_) / 2.0) * DiameterCorrectionParameter;
        sigma_m2 = M_PI * std::pow(correctedMeanDiameter, 2.0);
    }
    else{
        sigma_m2 = M_PI * std::pow(((ion.getDiameter() + collisionGasDiameter_m_) / 2.0), 2.0);
    }

    // Calculate the mean free path (MFP) from current ion velocity:

    // a static ion leads in static gas leads to a relative velocity of zero, which leads
    // to undefined behavior due to division by zero later.
    // The whole process converges to the MFP and collision probability of a static ion, thus
    // it is possible to assume a small velocity (1 nm/s) for the static ions to get rid of undefined behavior
    if (vRelIonMeanBackRest < 1e-9) {
        vRelIonMeanBackRest = 1e-9;
    }

    // Calculate the mean gas speed (m/s)
    double temperature_K = temperatureFunction_(pLocation);
    double vMeanGas = std::sqrt(
            8.0 * Core::K_BOLTZMANN * temperature_K / M_PI / (collisionGasMass_Amu_ * Core::AMU_TO_KG));

    // Calculate the median gas speed (m/s)
    double vMedianGas = std::sqrt(2.0 * Core::K_BOLTZMANN * temperature_K / (collisionGasMass_Amu_ * Core::AMU_TO_KG));

    // Compute the mean relative speed (m/s) between ion and gas.
    double s = vRelIonMeanBackRest / vMedianGas;
    double cMeanRel = vMeanGas * (
            (s + 1.0 / (2.0 * s)) * 0.5 * PI_SQRT * std::erf(s) + 0.5 * std::exp(-s * s));

    // Compute mean-free-path (m)
    double effectiveMFP_m = Core::K_BOLTZMANN * temperature_K *
                            (vRelIonMeanBackRest / cMeanRel) / (localPressure_Pa * sigma_m2);

    // Compute probability of collision in the current time-step.
    double collisionProb = 1.0 - std::exp(-vRelIonMeanBackRest * dt / effectiveMFP_m);

    // FIXME: The time step length dt is unrestricted
    // Possible mitigation: Throw warning / exception if collision probability becomes too high

    // Decide if a collision actually happens:
    if (rndSource->uniformRealRndValue() > collisionProb) {
        return; // no collision takes place
    }

    // Now we know that a collision happens: Perform the collision

     // Calculate the standard deviation of the one dimensional velocity distribution of the
    // background gas particles. Std. dev. in one dimension is given from Maxwell-Boltzmann
    // as sqrt(kT / particle mass).
    double vrStdevGas = std::sqrt(Core::K_BOLTZMANN * temperature_K / (collisionGasMass_Amu_ * Core::AMU_TO_KG));

    // Compute the velocity vector of the background gas particle colliding with the ion.

    Core::Vector vGasParticle;
    if (maxwellianApproximation_) {
        // Fast, approximate option: Calculate the gas particle velocity with a simple
        // Maxwell-Boltzmann distribution
        vGasParticle.x(rndSource->normalRealRndValue() * vrStdevGas);
        vGasParticle.y(rndSource->normalRealRndValue() * vrStdevGas);
        vGasParticle.z(rndSource->normalRealRndValue() * vrStdevGas);
    } else {
        // More correct but slower option:

        // A rejection method is used to account for the relative velocities between the ion and
        // the neutral background gas particles.
        // The probability of a gas particle with v_gas to hit the ion is
        // given by p(v_gas) = |v_gas - v_ion| f(v_gas) with f(v_gas) as the velocity
        // (Maxwell-Boltzmann) distribution of the background gas particles.
        double vGasParticleMagnitude;

        // vGasParticleUpperScale reasonable upper-bound for the length of the relative gas particle
        // velocity "vGasParticleMagnitude".
        // Here three standard deviations of the mean three dimensional background particle velocity is used
        double vGasParticleUpperScale = vRelIonMeanBackRest + vrStdevGas * SQRT3_3;
        do {
            vGasParticle.x(rndSource->normalRealRndValue() * vrStdevGas);
            vGasParticle.y(rndSource->normalRealRndValue() * vrStdevGas);
            vGasParticle.z(rndSource->normalRealRndValue() * vrStdevGas);
            vGasParticleMagnitude = (vGasParticle - vFrameMeanBackRest).magnitude();
        } while (rndSource->uniformRealRndValue() >= (vGasParticleMagnitude / vGasParticleUpperScale));
    }

    Core::Vector vIonPreCollision = ion.getVelocity();

    // Transform the frame of reference in a center of mass frame.
    double ionMass_kg = ion.getMass();
    double massDivisor = 1 / (ionMass_kg + collisionGasMass_kg_);
    Core::Vector cmfPreCollision;
    cmfPreCollision.x((ionMass_kg * vIonPreCollision.x() + collisionGasMass_kg_ * vGasParticle.x()) * massDivisor);
    cmfPreCollision.y((ionMass_kg * vIonPreCollision.y() + collisionGasMass_kg_ * vGasParticle.y()) * massDivisor);
    cmfPreCollision.z((ionMass_kg * vIonPreCollision.z() + collisionGasMass_kg_ * vGasParticle.z()) * massDivisor);

    // Get the alpha value
    double vss_collision_alpha = ion.getFloatAttribute(VSS_ALPHA);


    // Calculate the postcollision energy through E = 1/2 mv^2, with the velocity being the reference frame of
    // the colliding background gas particle at rest.
    Core::Vector vFrameCollidingBackRest = vFrameMeanBackRest - vGasParticle;
    // double reducedMass = ionMass_kg * collisionGasMass_kg_ / (ionMass_kg + collisionGasMass_kg_);
    // cannot differ from the precollsion energy (law of conservation of energy). Thus, the postcollision
    // energy can be calculated using the precollision viscosity vector.
    // double kineticEnergy = 0.5 * reducedMass * (vFrameMeanBackRest.magnitude()*vFrameMeanBackRest.magnitude());

    // Determine angle of the collision plane round the collision axis. All collision
    // planes are equally probable, since there is no preferential direction.
    double impactTheta = PI_2 * rndSource->uniformRealRndValue();

    // Calculate the resulting vectors
    // Core::Vector postCollisionVectorBackRest;

    // calculate cosX and sinX as impactAngle
    double cosX;
    double sinX;
    Core::Vector cmfPostCollision;

    // If-Else-Clauses to check if collision is approximately Hard Sphere.
    if (std::abs(1.0 - (1 / vss_collision_alpha)) < 0.001) {

        // Calculate the scattering angle, while assuming approximately Hard Sphere Collision. Thus, no dependency
        // on the alpha soft sphere scattering value.
        cosX = 2 * rndSource->uniformRealRndValue() - 1;

        double vFrameCollidingBackRestMagnitudeNeg = vFrameCollidingBackRest.magnitude() * (-1);

        cmfPostCollision.x(vFrameCollidingBackRestMagnitudeNeg * cosX);
        cmfPostCollision.y(vFrameCollidingBackRestMagnitudeNeg * std::cos(impactTheta));
        cmfPostCollision.z(vFrameCollidingBackRestMagnitudeNeg * std::sin(impactTheta));
    }
    else {
        cosX = 2 * std::pow(rndSource->uniformRealRndValue(), (1 / vss_collision_alpha)) - 1;
        sinX = std::sqrt(1 - cosX * cosX);

        double d = std::sqrt(cmfPreCollision.x() * cmfPreCollision.x() + cmfPreCollision.y() * cmfPreCollision.y());
        // Modify ion velocity in the normal direction due to elastic collision
        if (d < 1.0e-6){
            cmfPostCollision.x(cmfPreCollision.x() * cosX + sinX * d * std::sin(impactTheta));
            cmfPostCollision.y(cmfPreCollision.y() * cosX + (sinX * std::cos(impactTheta) * cmfPreCollision.z() -
                                                             cmfPreCollision.x() * cmfPreCollision.y() * sinX *
                                                             std::sin(impactTheta) / d));
            cmfPostCollision.z(cmfPreCollision.z() * cosX + (sinX * std::cos(impactTheta) * cmfPreCollision.y() -
                                                             cmfPreCollision.x() * cmfPreCollision.z() * sinX *
                                                             std::sin(impactTheta) / d));
        } else {
            cmfPostCollision.x(cmfPreCollision.x() * cosX);
            cmfPostCollision.y(cmfPreCollision.x() * std::cos(impactTheta));
            cmfPostCollision.z(cmfPreCollision.x() * std::sin(impactTheta));
        }
    }

    // Rotate reference frame back to the original reference frame
    Core::Vector vIonPostCollision;
    vIonPostCollision.x(cmfPreCollision.x() + (collisionGasMass_kg_ * massDivisor) * cmfPostCollision.x());
    vIonPostCollision.y(cmfPreCollision.y() + (collisionGasMass_kg_ * massDivisor) * cmfPostCollision.y());
    vIonPostCollision.z(cmfPreCollision.z() + (collisionGasMass_kg_ * massDivisor) * cmfPostCollision.z());

    // Set resulting velocity vector of the ion after the elastic collision
    ion.setVelocity(vIonPostCollision);


    // After the collision is finished:
    // Handle additional collision actions (e.g. collision based reactions etc.):

    // vFrameCollidingBackRest is the relative collision velocity
    // and thus the velocity for the reaction kinetic analysis
    if (afterCollisionActionFunction_ != nullptr) {
        // calculate collision energy in a center of mass frame
        // K_rel = 1/2 \mu v_rel^2$
        double m1 = ion.getMass();
        double m2 = collisionGasMass_kg_;
        double reducedMass = (m1 * m2) / (m1 + m2);
        double vRelativeMagnitude = vFrameCollidingBackRest.magnitude();
        double KEcollision = 0.5 * reducedMass * vRelativeMagnitude * vRelativeMagnitude;
        RS::CollisionConditions collisionConditions = {.totalCollisionEnergy = KEcollision};
        afterCollisionActionFunction_(collisionConditions, ion);
    }
}


void CollisionModel::SoftSphereModel::modifyPosition(Core::Vector& /*position*/, Core::Particle& /*ion*/, double /*dt*/) {}