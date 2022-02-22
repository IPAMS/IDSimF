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

#include "CollisionModel_StatisticalDiffusion.hpp"
#include "Core_randomGenerators.hpp"

/**
 * Constructs a statistical diffusion model with static background gas
 * (no velocity, static pressure and static temperature)
 *
 * @param staticPressure static background gas pressure (in Pa)
 * @param staticTemperature static background temperature (in K)
 * @param collisionGasMassAmu mass of the background gas particles (in amu)
 * @param collisionGasDiameterM effective collison diameter of the background gas particles (in m)
 * @param cs normalized diffusive collision statistics
 */
CollisionModel::StatisticalDiffusionModel::StatisticalDiffusionModel(
        double staticPressure,
        double staticTemperature,
        double collisionGasMassAmu,
        double collisionGasDiameterM,
        CollisionModel::CollisionStatistics cs):
        StatisticalDiffusionModel(
                staticPressure,
                staticTemperature,
                Core::Vector(0.0,0.0,0.0),
                collisionGasMassAmu,
                collisionGasDiameterM,
                std::move(cs)) {}

/**
 * Constructs a statistical diffusion model with uniformly flowing background gas
 * (static velocity vector, static pressure and static temperature)
 *
 * @param staticPressure static background gas pressure (in Pa)
 * @param staticTemperature static background temperature (in K)
 * @param staticGasVelocity static gas velocity (m/s)
 * @param collisionGasMassAmu mass of the background gas particles (in amu)
 * @param collisionGasDiameterM effective collison diameter of the background gas particles (in m)
 * @param cs normalized diffusive collision statistics
 */
CollisionModel::StatisticalDiffusionModel::StatisticalDiffusionModel(
        double staticPressure,
        double staticTemperature,
        Core::Vector staticGasVelocity,
        double collisionGasMassAmu,
        double collisionGasDiameterM,
        CollisionStatistics cs):
            StatisticalDiffusionModel(
                getConstantDoubleFunction(staticPressure),
                getConstantDoubleFunction(staticTemperature),
                getConstantVectorFunction(staticGasVelocity),
                collisionGasMassAmu,
                collisionGasDiameterM,
                std::move(cs))
{}

/**
 * Constructs a statistical diffusion model with variable background gas
 * (pressure, temperature and velocity given by arbitrary functions)
 *
 * @param pressureFunction Function to map spatial location to pressure (in Pa)
 * @param temperatureFunction Function to map spatial location to temperature (in K)
 * @param velocityFunction Function to map spatial location to velocity vector (in m/s)
 * @param collisionGasMassAmu mass of the background gas particles (in amu)
 * @param collisionGasDiameterM effective collison diameter of the background gas particles (in m)
 * @param cs normalized diffusive collision statistics
 */
CollisionModel::StatisticalDiffusionModel::StatisticalDiffusionModel(
        std::function<double(const Core::Vector &location)> pressureFunction,
        std::function<double(const Core::Vector &location)> temperatureFunction,
        std::function<Core::Vector(const Core::Vector &location)>velocityFunction,
        double collisionGasMassAmu,
        double collisionGasDiameterM,
        CollisionModel::CollisionStatistics cs):
            pressureFunction_(std::move(pressureFunction)),
            temperatureFunction_(std::move(temperatureFunction)),
            velocityFunction_(std::move(velocityFunction)),
            cs_(std::move(cs)),
            collisionGasMass_amu_(collisionGasMassAmu),
            collisionGasDiameter_nm_(collisionGasDiameterM * 1.0e9)
{
    // Init collision statistics:
    icdfs_ = cs_.getICDFs();
}

/**
 * Calculates a standard condition (normalized) random walk distance of a particle induced by diffusion.
 *
 * The actual collision number the diffusive random walk is representing is defined by the CollisionStatistics
 * object set for this collision model
 *
 * @param logParticleMassRatio decadic log of the mass ratio between particle and background gas particles
 * @return A normalized random walk distance
 */
double CollisionModel::StatisticalDiffusionModel::randomWalkDistance_(double logParticleMassRatio) const {
    // Find from the logarithmic mass ratio the indices of two ICFSs
    // in the collison statistics to interpolate between them.
    size_t upperDistIndex = cs_.findUpperDistIndex(logParticleMassRatio);

    // Get references to actual ICDFs (to prevent copying of ICFDs)
    const std::vector<double> &icdfLower = icdfs_.at(upperDistIndex);
    const std::vector<double> &icdfUpper = icdfs_.at(upperDistIndex+1);

    // Select actual indices in the discrete ICDFs according to an uniformly distributed random percentile
    double percentile = Core::globalRandomGeneratorPool->getThreadRandomSource()->uniformRealRndValue() * (cs_.getNDistPoints() - 2);
    size_t iLower = size_t(floor(percentile));
    size_t iUpper= iLower + 1;

    // Calculate random walk distance within the two ICDFs
    double weight = fmod(percentile,1.0);
    double rwdLower = (icdfLower[iUpper] - icdfLower[iLower]) * weight + icdfLower[iLower];
    double rwdUpper = (icdfUpper[iUpper] - icdfUpper[iLower]) * weight + icdfUpper[iLower];

    // Calculate a log interpolation between the ICDFs according to the logarithmic mass ratio (MR) of the particle
    rwdLower = log10(rwdLower);
    rwdUpper = log10(rwdUpper);

    // Calculate distance between log(particle MR) and log(MR) of lower ICDF,
    // normalize with distance between lower and upper ICDF mass ratio, then calculate final RWD
    double interIcdfWeight = (logParticleMassRatio - cs_.getLogMassRatio(upperDistIndex)) / cs_.getLogMassRatioDistance(upperDistIndex);
    double rwdFinal = pow(10,(rwdUpper - rwdLower) * interIcdfWeight + rwdLower);

    return rwdFinal;
}

/**
 * Calculates and sets the average velocity, the default stokes damping
 * and the average mean free path for an ion at standard conditions (STP)
 *
 * @param ion The ion to calculate the parameters for
 */
void CollisionModel::StatisticalDiffusionModel::setSTPParameters(BTree::Particle &ion) const{

    // FIXME: Cache the STP parameters instead of always recalculate them
    // (all parameters used here for the STP calculation are set / calculatable from the chemical species)
    // (create a cache structure with the stp values from the species as keys)

    double ionMass_kg = ion.getMass();
    double ionDiameter = ion.getDiameter();
    double ionDiameter_nm = ionDiameter * 1e9;

    // Calculate stokes damping factor (in s^-1):
    double damping =    Core::ELEMENTARY_CHARGE
                        / ion.getMobility()
                        / ionMass_kg;

    // Calculate a normalized mean speed from Maxwell Boltzmann distribution
    // (mean speed of a "particle" with unit mass )
    double vNormalized = sqrt(8 * Core::K_BOLTZMANN * STP_TEMP / M_PI);

    // Calculate thermal mean velocity of the ion (m/s)
    double vThermalIon = vNormalized * sqrt(1.0 / ionMass_kg);

    // Calculate thermal mean velocity of background gas particles (m/s)
    double vThermalGas = vNormalized * sqrt(1 / collisionGasMass_amu_ / Core::AMU_TO_KG);

    // Calculate collision frequency (collisions/s)
    double collisionFrequency = M2_PER_NM2 * STP_PARTICLE_DENSITY * M_PI *
                ((sqrt(2)-1.0/4.0) * pow(((collisionGasDiameter_nm_+ ionDiameter_nm)/2),2.0) * vThermalIon +
                 (1.0/4.0) * pow(ionDiameter_nm,2.0) * vThermalGas);

    // Calculate mean free path (mfp) of the ion at STP (in m)
    double mfp = vThermalIon / collisionFrequency;

    ion.setMeanThermalVelocitySTP(vThermalIon);
    ion.setMeanFreePathSTP(mfp);
    ion.getAuxCollisionParams()[index_ionSTPdamping] = damping;
}

/**
 * Updates internal model parameters for an ion which are dependent on ion position / timestep
 */
void CollisionModel::StatisticalDiffusionModel::updateModelParameters(BTree::Particle &ion) const {

    Core::Vector particle_location = ion.getLocation();

    // Get temperature and pressure at the ion location
    double localTemperature_K = temperatureFunction_(particle_location);
    double localPressure_pa  = pressureFunction_(particle_location);

    // Calculate temperature and pressure/temperature ratios for the ion
    double tRatio = localTemperature_K / STP_TEMP;
    double ptRatio = tRatio * (STP_PRESSURE / localPressure_pa);

    // Set parameters for the ion (ratios for the ion are used in other methods)
    ion.getAuxCollisionParams()[index_tRatio] = tRatio;
    ion.getAuxCollisionParams()[index_ptRatio] = ptRatio;
}

/**
 * Inits model parameters which are not dependent on ion position / timestep for an ion
 */
void CollisionModel::StatisticalDiffusionModel::initializeModelParameters(BTree::Particle &ion) const {
    this->setSTPParameters(ion);
}

/**
 * Modifies the acceleration due to the background gas interaction
 *
 * @param acceleration The current acceleration of the ion which is modified (the acceleration
 * stored in the given particle object is usually later overwritten in the integration mechanism)
 * @param ion An ion the acceleration is calculated for
 * @param dt The length of the current time step
 */
void CollisionModel::StatisticalDiffusionModel::modifyAcceleration(Core::Vector& acceleration, BTree::Particle& ion,
                                                                   double dt) {

    Core::Vector gasVelocity = velocityFunction_(ion.getLocation());
    double ionLocalDamping = ion.getAuxCollisionParams()[index_ionSTPdamping]/
            ion.getAuxCollisionParams()[index_ptRatio];

    double dampingTimeConstant = ionLocalDamping * dt;
    double dampingFactor = (1 - exp(-dampingTimeConstant)) / dampingTimeConstant;

    // Calculate new acceleration (a = dv / dt)
    Core::Vector newAcceleration =
            (acceleration - (ion.getVelocity() - gasVelocity) * ionLocalDamping )
            * dampingFactor;

    acceleration = newAcceleration;
}

/**
 * The ion velocity is not modified by the SDS model
 */
void CollisionModel::StatisticalDiffusionModel::modifyVelocity(BTree::Particle& /*ion*/, double /*dt*/) {}

/**
 * Modifies the position of a particle according to the background gas interaction modeled with the SDS approach
 *
 * @param position The current position of the particle which is modified (the position
 * stored in the particle object is usually later overwritten in the integration mechanism)
 * @param ion A ion the position is calculated for
 * @param dt The length of the current time step
  */
void CollisionModel::StatisticalDiffusionModel::modifyPosition(
        Core::Vector& position, BTree::Particle &ion, double dt) {

    Core::Vector oldPosition = position;
    double ionMass_amu = ion.getMass() / Core::AMU_TO_KG;
    double t_ratio = ion.getAuxCollisionParams()[index_tRatio];
    double pt_ratio = ion.getAuxCollisionParams()[index_ptRatio];

    // Calculate mass ratio between ions and background gas
    // and current mean free path and thermal velocity for the ion:
    double massRatioLog = log10(ionMass_amu / collisionGasMass_amu_);
    double ionMfpLocal = ion.getMeanFreePathSTP() * pt_ratio;
    double ionThermalVelocityLocal = ion.getMeanThermalVelocitySTP() * sqrt(t_ratio);

    // Calculate normalized RWD and the average number of collisions the ion experiences:
    double randomWalkDistanceNormalized = randomWalkDistance_(massRatioLog);
    double nCollisions=
            ionThermalVelocityLocal  // (m/s)
            / ionMfpLocal            // (1/m)
            * dt;                    // (s)

    // Calculate the actual RWD for the conditions of the ion:
    double randomJumpDistance =
            sqrt(nCollisions / cs_.getNCollisions())  // Scale between actual coll. and coll. in coll. statistics
            * randomWalkDistanceNormalized            // Normalized RWD with mean free path = 1
            * ionMfpLocal;                            // Local MFP of the ion

    // Apply random jump:
    Core::Vector randomJump = sphereRand(randomJumpDistance);
    position = oldPosition + randomJump;
}