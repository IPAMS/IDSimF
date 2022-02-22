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
 PSim_velocityIntegrator.hpp

 Trajectory integrator which takes just a simple velocity function for the particles
 (no space charge or collision simulation)

 ****************************/

// FIXME: Velocity integrator is currently not able to have Time Of Birth Distributions / Bear Particles

#ifndef BTree_velocityIntegrator_hpp
#define BTree_velocityIntegrator_hpp

#include "PSim_abstractTimeIntegrator.hpp"
#include "Core_vector.hpp"
#include <cstdio>
#include <vector>
#include <functional>
//forward declare own classes:
namespace BTree{
    class Particle;
}

namespace Integration{

    /**
     * Ion trajectory integration with a very simple integrator which moves particles in every timestep according to
     * their local velocity.
     *
     * The velocity calculation and what is exported in every time step is passed to this trajectory integrator
     * externally by functions. Thus, the integration scheme can be applied to various simulation problems.
     *
     * No space charge or background collision interaction is considered with this integrator.
     */
    class VelocityIntegrator: public AbstractTimeIntegrator{

    public:

        /**
        * type definition for velocity calculation functions
        */
        typedef std::function
                <Core::Vector (BTree::Particle* particle,
                                std::size_t particleIndex,
                                double time,
                                unsigned int timestep)>
                velocityFctType;

        /**
         * type definition for functions exporting data in every timestep
         */
        typedef std::function
                <void (std::vector<BTree::Particle*>& particles,
                       double time,
                       unsigned int timestep,
                       bool lastTimestep)>
                timestepWriteFctType;

        /**
         * type definition for functions defining "other actions", which are additional arbitrary actions performed
         * in every time step of the integration
         */
        typedef std::function
                <void (Core::Vector& newPartPos,
                       BTree::Particle* particle,
                       std::size_t particleIndex,
                       double time,
                       unsigned int timestep)>
                otherActionsFctType;

        VelocityIntegrator(
                std::vector<BTree::Particle*> particles,
                velocityFctType velocityFunction,
                timestepWriteFctType timestepWriteFunction = nullptr,
                otherActionsFctType otherActionsFunction = nullptr
        );

        void addParticle(BTree::Particle* particle) override;
        void run(unsigned int nTimesteps, double dt) override;
        void runSingleStep(double dt) override;
        void finalizeSimulation() override;

    private:
        void addParticle_(BTree::Particle* particle);

        velocityFctType velocityFunction_ = nullptr; ///< function to calculate particle acceleration
        timestepWriteFctType timestepWriteFunction_ = nullptr; ///< function to define what is exported in every time step
        otherActionsFctType otherActionsFunction_ = nullptr; ///< function with other actions done in every time step
    };
}
#endif /* BTree_velocityIntegrator_hpp */