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
 Integration_fmm3dIntegrator.hpp

 Fast multipole method (FMM) based ion trajectory integration
 (simple verlet trajectory integrator) including space charge

 ****************************/

#ifndef IDSIMF_INTEGRATION_FMM3DINTEGRATOR_HPP
#define IDSIMF_INTEGRATION_FMM3DINTEGRATOR_HPP

#include "Core_particle.hpp"
#include "Core_vector.hpp"
#include "Integration_abstractTimeIntegrator.hpp"
#include "CollisionModel_AbstractCollisionModel.hpp"
#include "FMM3D_fmmSolver.hpp"
#include <vector>
namespace Integration{

    //std::function<Core::Vector(Core::Particle* particle, int particleIndex, Core::Tree& tree, double time, int timestep)> accelerationFctType;

    class FMMVerletIntegrator : public AbstractTimeIntegrator {
    public:

        FMMVerletIntegrator(
            const std::vector<Core::Particle*>& particles,
            accelerationFctType accelerationFunction,
            timestepWriteFctType timestepWriteFunction = nullptr,
            otherActionsFctType otherActionsFunction = nullptr,
            AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction = nullptr,
            CollisionModel::AbstractCollisionModel* collisionModel = nullptr
        );

        FMMVerletIntegrator(
            accelerationFctType accelerationFunction,
            timestepWriteFctType timestepWriteFunction = nullptr,
            otherActionsFctType otherActionsFunction = nullptr,
            AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction = nullptr,
            CollisionModel::AbstractCollisionModel* collisionModel = nullptr
        );

        void addParticle(Core::Particle* particle) override;
        void run(unsigned int nTimesteps, double dt) override;
        void runSingleStep(double dt) override;
        void finalizeSimulation() override;

    private:
        FMM3D::FMMSolver solver_;

        CollisionModel::AbstractCollisionModel* collisionModel_ = nullptr; ///< the gas collision model to perform while integrating

        accelerationFctType accelerationFunction_ = nullptr;   ///< function to calculate particle acceleration
        timestepWriteFctType timestepWriteFunction_ = nullptr; ///< function to export / write time step results
        otherActionsFctType otherActionsFunction_ = nullptr;   ///< function for arbitrary other actions in the simulation

        std::vector<Core::Vector>  a_t_;     ///< last time step acceleration for particles
        std::vector<Core::Vector>  a_tdt_;   ///< new acceleration for particles

    };
}


#endif //IDSIMF_INTEGRATION_FMM3DINTEGRATOR_HPP
