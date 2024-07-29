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
 BTree_fullSumRK4Integrator.hpp
 
 Parallel version of a Runge Kutta 4 Integrator with naive full sum space charge calculation
 (mostly for benchmarking and testing purposes)

 ****************************/

#ifndef BTree_fullSumRK4Integrator_hpp
#define BTree_fullSumRK4Integrator_hpp

#include "Core_particle.hpp"
#include "Core_vector.hpp"
#include "SC_fullSumSolver.hpp"
#include "Integration_abstractTimeIntegrator.hpp"
#include "CollisionModel_AbstractCollisionModel.hpp"
#include <vector>

namespace Integration{

    class FullSumRK4Integrator: public AbstractTimeIntegrator {

        public:

            FullSumRK4Integrator(
                    const std::vector<Core::Particle*>& particles,
                    accelerationFctType accelerationFunction,
                    accelerationFctSpaceChargeType spaceChargeAccelerationFunction,
                    postTimestepFctType timestepWriteFunction = nullptr,
                    otherActionsFctType otherActionsFunction = nullptr,
                    AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction = nullptr,
                    CollisionModel::AbstractCollisionModel* collisionModel = nullptr
            );

            FullSumRK4Integrator(
                    accelerationFctType accelerationFunction,
                    accelerationFctSpaceChargeType spaceChargeAccelerationFunction,
                    postTimestepFctType timestepWriteFunction = nullptr,
                    otherActionsFctType otherActionsFunction = nullptr,
                    AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction = nullptr,
                    CollisionModel::AbstractCollisionModel* collisionModel = nullptr
            );

            void addParticle(Core::Particle* particle) override;
            void run(unsigned int nTimesteps, double dt) override;
            void runSingleStep(double dt) override;
            void finalizeSimulation() override;

    private:

        CollisionModel::AbstractCollisionModel* collisionModel_ = nullptr; ///< the gas collision model to perform while integrating

        accelerationFctType accelerationFunction_ = nullptr;   ///< function to calculate particle acceleration (without space charge)
        accelerationFctSpaceChargeType spaceChargeAccelerationFunction_ = nullptr;   ///< function to calculate particle acceleration from space charge
        postTimestepFctType postTimestepFunction_ = nullptr; ///< function to export / write time step results
        otherActionsFctType otherActionsFunction_ = nullptr;   ///< function for arbitrary other actions in the simulation

        Core::Vector evaluateAccelerationFunction_(Core::Particle* particle,
                                                   Core::Vector position,
                                                   Core::Vector velocity,
                                                   Core::Vector spaceChargeAcceleration,
                                                   double time,
                                                   unsigned int timestep,
                                                   double dt);

        //internal variables for actual calculations:

        std::vector<Core::Vector>  newPos_;  ///< new position (after time step) for particles
        std::vector<Core::Vector>  spaceChargeAcceleration_;  ///< space charge acceleration from

        SpaceCharge::FullSumSolver fullSumSolver_ = SpaceCharge::FullSumSolver();  ///< full sum space charge solver
        void bearParticles_(double time);
    };
}


#endif /* BTree_fullSumRK4Integrator_hpp */
