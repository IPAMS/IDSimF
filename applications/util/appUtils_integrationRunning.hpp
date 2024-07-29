/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2024 - Physical and Theoretical Chemistry /
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
 appUtils_integrationRunning.hpp

 Helping / Convenience Methods for simple simulation / integration running

 ****************************/
#ifndef IDSIMF_APPUTILS_INTEGRATIONRUNNING_HPP
#define IDSIMF_APPUTILS_INTEGRATIONRUNNING_HPP

#include "Core_particle.hpp"
#include "Integration_generic.hpp"
#include "Integration_abstractTimeIntegrator.hpp"
#include "CollisionModel_AbstractCollisionModel.hpp"
#include "appUtils_simulationConfiguration.hpp"

namespace AppUtils{

    bool integratorModeInVector(const std::vector<AppUtils::IntegratorMode>& vec, AppUtils::IntegratorMode mode);

    void runTrajectoryIntegration(
            AppUtils::simConf_ptr simConf,
            unsigned int timeSteps,
            double dt,
            const std::vector<Core::Particle*>& particles,
            Integration::accelerationFctSingleStepType singleStepAccelerationFunction,
            Integration::postTimestepFctType postTimestepFunction = nullptr,
            Integration::otherActionsFctType otherActionsFunction = nullptr,
            Integration::AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction = nullptr,
            CollisionModel::AbstractCollisionModel* collisionModel = nullptr);

    void runTrajectoryIntegration(
            AppUtils::simConf_ptr simConf,
            unsigned int timeSteps,
            double dt,
            const std::vector<Core::Particle*>& particles,
            Integration::accelerationFctType multiStepAccelerationFunction,
            Integration::accelerationFctSpaceChargeType multiStepSpaceChargeAccelerationFunction,
            Integration::postTimestepFctType postTimestepFunction = nullptr,
            Integration::otherActionsFctType otherActionsFunction = nullptr,
            Integration::AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction = nullptr,
            CollisionModel::AbstractCollisionModel* collisionModel = nullptr);

    void runTrajectoryIntegration(
            AppUtils::simConf_ptr simConf,
            unsigned int timeSteps,
            double dt,
            const std::vector<Core::Particle*>& particles,
            Integration::accelerationFctSingleStepType singleStepAccelerationFunction,
            Integration::accelerationFctType multiStepAccelerationFunction,
            Integration::accelerationFctSpaceChargeType multiStepSpaceChargeAccelerationFunction,
            Integration::postTimestepFctType postTimestepFunction = nullptr,
            Integration::otherActionsFctType otherActionsFunction = nullptr,
            Integration::AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction = nullptr,
            CollisionModel::AbstractCollisionModel* collisionModel = nullptr);
}

#endif //IDSIMF_APPUTILS_INTEGRATIONRUNNING_HPP
