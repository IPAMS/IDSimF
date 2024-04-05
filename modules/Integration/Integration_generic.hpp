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

 ------------
 Integration_generic.hpp

 Generic interfaces / type definitions for time integration

 ****************************/
#ifndef IDSIMF_INTEGRATION_GENERIC_HPP
#define IDSIMF_INTEGRATION_GENERIC_HPP

#include "SC_generic.hpp"
#include <vector>
#include <functional>

namespace Integration{

    // Forward abstract time integrator
    class AbstractTimeIntegrator;

    /**
     * type definition for acceleration calculation functions (single step velocity verlet typed acceleration)
     */
    typedef std::function
        <Core::Vector (
            Core::Particle* particle,
            std::size_t particleIndex,
            SpaceCharge::FieldCalculator &forceCalculator,
            double time,
            unsigned int timestep)>
    accelerationFctSingleStepType;

    /**
     * Acceleration function for arbitrary position, velocity, time and mass
     * (typically for pure forces without space charge in multistep methods)
     */
    typedef std::function
            <Core::Vector (
                    Core::Particle* particle,
                    Core::Vector position,
                    Core::Vector velocity,
                    double time,
                    unsigned int timestep)>
    accelerationFctType;

    /**
     * Acceleration of a particle from space charge
     * (typically for space charge contribution in multistep integrator methods)
     */
    typedef std::function
            <Core::Vector (
                    Core::Particle* particle,
                    std::size_t particleIndex,
                    SpaceCharge::FieldCalculator &forceCalculator,
                    double time,
                    unsigned int timestep)>
    accelerationFctSpaceChargeType;



    /**
     * type definition for functions doing things after every timestep (mostly exporting data in every timestep or
     * stopping the integration)
     */
    typedef std::function
        <void (
            AbstractTimeIntegrator* integrator,
            std::vector<Core::Particle*>& particles,
            double time,
            unsigned int timestep,
            bool lastTimestep)>
    postTimestepFctType;

    /**
     * type definition for functions defining "other actions", which are additional arbitrary actions performed
     * in every time step of the integration
     */
    typedef std::function
        <void (
            Core::Vector& newPartPos,
            Core::Particle* particle,
            std::size_t particleIndex,
            double time,
            unsigned int timestep)>
    otherActionsFctType;

}

#endif //IDSIMF_INTEGRATION_GENERIC_HPP
