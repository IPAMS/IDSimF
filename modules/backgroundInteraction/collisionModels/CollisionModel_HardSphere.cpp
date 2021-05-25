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

#include "CollisionModel_HardSphere.hpp"
#include "Core_math.hpp"
#include "Core_randomGenerators.hpp"

/**
 * Constructor for static background gas pressure and temperature
 *
 * if maxwellianApproximation is true, the colliding background gas particle velocity is drawn
 * from a Maxwell Boltzmann distribution instead of the more correct relative velocity distribution
 * between charged particle and background gas
 */
CollisionModel::HardSphereModel::HardSphereModel(
        double staticPressure,
        double staticTemperature,
        double collisionGasMassAmu,
        double collisionGasDiameterM,
        bool maxwellianApproximation)
        :
        HardSphereModel(
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
CollisionModel::HardSphereModel::HardSphereModel(
        double staticPressure,
        double staticTemperature,
        double collisionGasMassAmu,
        double collisionGasDiameterM,
        std::function<void(RS::CollisionConditions, BTree::Particle&)> afterCollisionFunction,
        bool maxwellianApproximation)
        :
        HardSphereModel(
                getConstantDoubleFunction(staticPressure),
                getConstantVectorFunction(Core::Vector(0.0, 0.0, 0.0)),
                staticTemperature,
                collisionGasMassAmu,
                collisionGasDiameterM,
                afterCollisionFunction,
                maxwellianApproximation) { }

/**
 * Constructor for location dependent pressure and velocity, given as spatially resolved
 * functions.
 */
CollisionModel::HardSphereModel::HardSphereModel(
        std::function<double(Core::Vector& location)> pressureFunction,
        std::function<Core::Vector(Core::Vector& location)> velocityFunction,
        double staticTemperature,
        double collisionGasMassAmu,
        double collisionGasDiameterM,
        bool maxwellianApproximation)
        :
        maxwellianApproximation_(maxwellianApproximation),
        temperature_K_(staticTemperature),
        collisionGasMass_Amu_(collisionGasMassAmu),
        collisionGasMass_kg_(collisionGasMassAmu*Core::AMU_TO_KG),
        collisionGasDiameter_m_(collisionGasDiameterM),
        pressureFunction_(std::move(pressureFunction)),
        velocityFunction_(std::move(velocityFunction)) { }

/**
 * Constructor for location dependent pressure and velocity, given as spatially resolved
 * functions and custom after collision function (most probably to model collision based reactions).
 */
CollisionModel::HardSphereModel::HardSphereModel(
        std::function<double(Core::Vector& location)> pressureFunction,
        std::function<Core::Vector(Core::Vector& location)> velocityFunction,
        double staticTemperature,
        double collisionGasMassAmu,
        double collisionGasDiameterM,
        std::function<void(RS::CollisionConditions, BTree::Particle&)> afterCollisionFunction,
        bool maxwellianApproximation)
        :
        maxwellianApproximation_(maxwellianApproximation),
        temperature_K_(staticTemperature),
        collisionGasMass_Amu_(collisionGasMassAmu),
        collisionGasMass_kg_(collisionGasMassAmu*Core::AMU_TO_KG),
        collisionGasDiameter_m_(collisionGasDiameterM),
        pressureFunction_(std::move(pressureFunction)),
        velocityFunction_(std::move(velocityFunction)),
        afterCollisionActionFunction_(std::move(afterCollisionFunction)) { }

void CollisionModel::HardSphereModel::updateModelParameters(BTree::Particle& /*ion*/) const {}

void CollisionModel::HardSphereModel::initializeModelParameters(BTree::Particle& /*ion*/) const {}

void CollisionModel::HardSphereModel::modifyAcceleration(Core::Vector& /*acceleration*/, BTree::Particle& /*ion*/,
                                                         double /*dt*/) {}

 /**
  * Modify the velocity of a charged particle "ion" with a potential random collision in time step "dt"
  * (only one or none collision can happen, not multiple. Thus, the probability of multiple colllision events
  * happening in "dt" has to be low)
  */
void CollisionModel::HardSphereModel::modifyVelocity(BTree::Particle &ion, double dt) {

    // Calculate collision cross section between particle and collision gas:
    //   TODO: It seems to be unnecessary to constantly recalculate this
    //   value, cache the calculated values somehow?
    double sigma_m2 = M_PI * std::pow( (ion.getDiameter() + collisionGasDiameter_m_)/2.0, 2.0);

    Core::Vector pLocation = ion.getLocation();
    double localPressure_Pa = pressureFunction_(pLocation);

    if (localPressure_Pa==0.0){
        return; //pressure 0 means no collision at all
    }

    // Transform the frame of reference in a frame where the mean background gas velocity is zero.
    Core::Vector vGasMean = velocityFunction_(pLocation);
    Core::Vector vFrameMeanBackRest = ion.getVelocity() - vGasMean;

    double vRelIonMeanBackRest = vFrameMeanBackRest.magnitude(); //relative ion relative to bulk gas velocity

    // Calculate the mean free path (MFP) from current ion velocity:

    // a static ion leads in static gas leads to a relative velocity of zero, which leads
    // to undefined behavior due to division by zero later.
    // The whole process converges to the MFP and collision probability of a static ion, thus
    // it is possible to assume a small velocity (1 nm/s) for the static ions to get rid of undefined behavior
    if (vRelIonMeanBackRest < 1e-9){
        vRelIonMeanBackRest = 1e-9;
    }

    // Calculate the mean gas speed (m/s)
    double vMeanGas = std::sqrt(8.0*Core::K_BOLTZMANN*temperature_K_/M_PI/(collisionGasMass_Amu_ * Core::AMU_TO_KG));

    // Calculate the median gas speed (m/s)
    double vMedianGas = std::sqrt(2.0*Core::K_BOLTZMANN*temperature_K_/(collisionGasMass_Amu_ * Core::AMU_TO_KG));

    // Compute the mean relative speed (m/s) between ion and gas.
    double s = vRelIonMeanBackRest / vMedianGas;
    double cMeanRel = vMeanGas * (
            (s + 1.0/(2.0*s)) * 0.5 * PI_SQRT * std::erf(s) + 0.5 * std::exp(-s*s) );

    // Compute mean-free-path (m)
    double effectiveMFP_m = Core::K_BOLTZMANN * temperature_K_ *
                            (vRelIonMeanBackRest / cMeanRel) / (localPressure_Pa * sigma_m2);

    // Compute probability of collision in the current time-step.
    double collisionProb = 1.0 - std::exp(-vRelIonMeanBackRest * dt / effectiveMFP_m);

    // FIXME: The time step length dt is unrestricted
    // Possible mitigation: Throw warning / exception if collision probability becomes too high

    // Decide if a collision actually happens:
    if (Core::globalRandomGenerator->uniformRealRndValue() > collisionProb){
        return; // no collision takes place
    }

    // Now we know that a collision happens: Perform the collision

    // Calculate the standard deviation of the one dimensional velocity distribution of the
    // background gas particles. Std. dev. in one dimension is given from Maxwell-Boltzmann
    // as sqrt(kT / particle mass).
    double  vrStdevGas = std::sqrt( Core::K_BOLTZMANN * temperature_K_ / (collisionGasMass_Amu_ * Core::AMU_TO_KG) );

    // Compute the velocity vector of the background gas particle colliding with the ion.

    Core::Vector vGasParticle;
    if (maxwellianApproximation_){
        // Fast, approximate option: Calculate the gas particle velocity with a simple
        // Maxwell-Boltzmann distribution
        vGasParticle.x(Core::globalRandomGenerator->normalRealRndValue() * vrStdevGas);
        vGasParticle.y(Core::globalRandomGenerator->normalRealRndValue() * vrStdevGas);
        vGasParticle.z(Core::globalRandomGenerator->normalRealRndValue() * vrStdevGas);
    }
    else {
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
            vGasParticle.x(Core::globalRandomGenerator->normalRealRndValue() * vrStdevGas);
            vGasParticle.y(Core::globalRandomGenerator->normalRealRndValue() * vrStdevGas);
            vGasParticle.z(Core::globalRandomGenerator->normalRealRndValue() * vrStdevGas);
            vGasParticleMagnitude = (vGasParticle - vFrameMeanBackRest).magnitude();
        }
        while (Core::globalRandomGenerator->uniformRealRndValue() >= (vGasParticleMagnitude / vGasParticleUpperScale));
    }




    // Define a new reference frame with the colliding background gas particle at rest
    // for the subsequent analysis
    Core::Vector vFrameCollidingBackRest = vFrameMeanBackRest - vGasParticle;

    // The collision between the two particles is handled in the 2d collision plane
    // (a plane in which the trajectories of the two collision remains during
    // the collision event)
    // The local coordinate system (u,v) in this plane is the collision axis u between the
    // colliding particles and the perpendicular axis v. The whole plane is rotated around
    // u by a rotation angle theta

    // Compute the impact offset of the collision.
    // impactOffset = 0 means a head-on collision while 1 means the particles just missed.
    // The impactOffset correspond to a ring of possible collision positions around the collision
    // axis. f(d) = 2*d with d in an interval [0,1] is the normalized probability density function
    // for a hit on a ring with radius d. This can be transformed with the fundamental transformation
    // law of probabilities ("Change of variables formula") to impact_offset = sqrt(Uni) with an
    // uniformly distributed random variable in the interval [0,1].
    double impactOffset = std::sqrt(Core::globalRandomGenerator->uniformRealRndValue());

    // Calculate the impact angle on the surface of spherical gas particle where the ion hits, which is
    // the vector from the gas particle center to the collision point on the surface
    double impactAngle = std::asin(impactOffset);

    // Determine angle of the collision plane round the collision axis. All collision
    // planes are equally probable, since there is no preferential direction.
    double impactTheta = PI_2*Core::globalRandomGenerator->uniformRealRndValue();

    // Compute spherical coordinates in current velocity reference frame.
    Core::Vector vFrameCollidingRest_sp = cartesianToPolar(vFrameCollidingBackRest);

    // (the resulting vector contains now polar coordinates in degrees)
    double vIon_sp = vFrameCollidingRest_sp.x();
    double azimuthIon = vFrameCollidingRest_sp.y();
    double elevationIon = vFrameCollidingRest_sp.z();

    // Velocity components of the ion relative to the collision axis
    // (connection axis between the particle centers in the moment of collision)
    double vIonNormal = vIon_sp * cos(impactAngle);   //normal velocity
    double vIonRadial = vIon_sp * sin(impactAngle);   //radial velocity

    // Modify ion velocity in the normal direction due to elastic collision
    // The force acts in the normal direction, which is also normal to the collision plane
    // (Tangential plane of the particle surfaces in the point of contact)
    double ionMass_kg = ion.getMass();
    double vIonNormalAfterCollision = (vIonNormal * (ionMass_kg - collisionGasMass_kg_))
                    / (ionMass_kg + collisionGasMass_kg_);

    // Rotate frame in a way that ion velocity is on the y axis, which also means
    // that the angle between the y axis and the resulting velocity vector is the
    // angle the particle was scattered by
    Core::Vector vFrameRot = elevationRotate(
            Core::Vector(vIonNormalAfterCollision, vIonRadial, 0),
            M_PI_2 - impactAngle);

    // Select the orientation of the plane the collision is taking place in by
    // rotating around the y axis with the angle theta
    vFrameRot = azimuthRotate(vFrameRot, impactTheta);

    // Rotate reference frame back to the original reference frame
    vFrameRot = elevationRotate(vFrameRot, -M_PI_2 + elevationIon);
    vFrameRot = azimuthRotate(vFrameRot, azimuthIon);

    // Translate reference frame back to original velocity rference frame
    // relative to the ion
    Core::Vector vFrameFinal = vFrameRot + vGasParticle + vGasMean;

    // Set resulting velocity vector of the ion after the elastic collision
    ion.setVelocity(vFrameFinal);

    // After the collision is finished:
    // Handle additional collision actions (e.g. collision based reactions etc.):

    // vFrameCollidingBackRest is the relative collision velocity
    // and thus the velocity for the reaction kinetic analysis
    if (afterCollisionActionFunction_ != nullptr){
        // calculate collision energy in a center of mass frame
        // K_rel = 1/2 \mu v_rel^2$
        double m1 = ion.getMass();
        double m2 = collisionGasMass_kg_;
        double reducedMass = (m1*m2)/(m1+m2);
        double vRelativeMagnitude = vFrameCollidingBackRest.magnitude();
        double KEcollision = 0.5 * reducedMass * vRelativeMagnitude * vRelativeMagnitude;
        RS::CollisionConditions collisionConditions = {.totalCollisionEnergy = KEcollision};
        afterCollisionActionFunction_(collisionConditions, ion);
    }
}

void CollisionModel::HardSphereModel::modifyPosition(Core::Vector& /*position*/, BTree::Particle& /*ion*/, double /*dt*/) {}