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
 appUtils_integrationRunning.cpp

 Description

 ****************************/
#include "appUtils_integrationRunning.hpp"
#include "appUtils_signalHandler.hpp"
#include "Integration_verletIntegrator.hpp"
#include "Integration_parallelVerletIntegrator.hpp"
#include "Integration_fullSumVerletIntegrator.hpp"
#include "Integration_parallelRK4Integrator.hpp"
#include "Integration_fullSumRK4Integrator.hpp"

#ifdef WITH_FMM_3d
#include "Integration_fmmIntegrator.hpp"
#include "FMM3D_fmmSolver.hpp"
#endif

#ifdef WITH_EXAFMMT
#include "Integration_fmmIntegrator.hpp"
#include "ExaFMMt_fmmSolver.hpp"
#endif

#include <stdexcept>

bool AppUtils::integratorModeInVector(const std::vector<AppUtils::IntegratorMode>& vec, AppUtils::IntegratorMode mode) {
    return std::find(vec.begin(), vec.end(), mode) != vec.end();
}

/**
 * Trajectory integrator runner for single step integrators only
 */
void AppUtils::runTrajectoryIntegration(
        AppUtils::simConf_ptr simConf,
        unsigned int timeSteps,
        double dt,
        const std::vector<Core::Particle*>& particles,
        Integration::accelerationFctSingleStepType singleStepAccelerationFunction,
        Integration::postTimestepFctType postTimestepFunction,
        Integration::otherActionsFctType otherActionsFunction,
        Integration::AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction,
        CollisionModel::AbstractCollisionModel* collisionModel) {

    IntegratorMode integratorMode = simConf->integratorMode();

    std::vector<IntegratorMode> singleStepIntegrators =
            {IntegratorMode::VERLET, IntegratorMode::PARALLEL_VERLET, IntegratorMode::FULL_SUM_VERLET,
             IntegratorMode::EXAFMM_VERLET, IntegratorMode::FMM3D_VERLET};

    if (integratorModeInVector(
            singleStepIntegrators,
            integratorMode)) {

        runTrajectoryIntegration(
                simConf, timeSteps, dt,
                particles,
                singleStepAccelerationFunction,
                nullptr, nullptr,
                postTimestepFunction, otherActionsFunction, ionStartMonitoringFunction, collisionModel);
    }
    else {
        throw (std::invalid_argument(
                "Wrong Integrator: A multistep integrator was selected for a simulation app which implements only single step integrators"));
    }
}


/**
 * Trajectory integrator runner for multi step integrators
 */
void AppUtils::runTrajectoryIntegration(
        AppUtils::simConf_ptr simConf,
        unsigned int timeSteps,
        double dt,
        const std::vector<Core::Particle*>& particles,
        Integration::accelerationFctType multiStepAccelerationFunction,
        Integration::accelerationFctSpaceChargeType multiStepSpaceChargeAccelerationFunction,
        Integration::postTimestepFctType postTimestepFunction,
        Integration::otherActionsFctType otherActionsFunction,
        Integration::AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction,
        CollisionModel::AbstractCollisionModel* collisionModel){

    IntegratorMode integratorMode = simConf->integratorMode();
    std::vector<IntegratorMode> multiStepIntegrators = {IntegratorMode::PARALLEL_RUNGE_KUTTA4, IntegratorMode::FULL_SUM_RUNGE_KUTTA4};
    if (AppUtils::integratorModeInVector(
            multiStepIntegrators, integratorMode)) {

        runTrajectoryIntegration(
                simConf, timeSteps, dt,
                particles,
                nullptr,
                multiStepAccelerationFunction, multiStepSpaceChargeAccelerationFunction,
                postTimestepFunction, otherActionsFunction, ionStartMonitoringFunction, collisionModel);
    }
    else {
        throw (std::invalid_argument(
                "Wrong Integrator: A single step integrator was selected for a simulation app which implements only multi step integrators"));
    }
}


/**
 * Trajectory integrator runner for single and multi step integrators
 */
void AppUtils::runTrajectoryIntegration(
        AppUtils::simConf_ptr simConf,
        unsigned int timeSteps,
        double dt,
        const std::vector<Core::Particle*>& particles,
        Integration::accelerationFctSingleStepType singleStepAccelerationFunction,
        Integration::accelerationFctType multiStepAccelerationFunction,
        Integration::accelerationFctSpaceChargeType multiStepSpaceChargeAccelerationFunction,
        Integration::postTimestepFctType postTimestepFunction,
        Integration::otherActionsFctType otherActionsFunction,
        Integration::AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction,
        CollisionModel::AbstractCollisionModel* collisionModel) {

    AppUtils::IntegratorMode integratorMode = simConf->integratorMode();

    if (integratorMode==AppUtils::VERLET) {
        Integration::VerletIntegrator verletIntegrator(
                particles,
                singleStepAccelerationFunction, postTimestepFunction, otherActionsFunction,
                ionStartMonitoringFunction,
                collisionModel);
        AppUtils::SignalHandler::setReceiver(verletIntegrator);
        verletIntegrator.run(timeSteps, dt);
    }
    else if (integratorMode==AppUtils::PARALLEL_VERLET) {
        Integration::ParallelVerletIntegrator verletIntegrator(
                particles,
                singleStepAccelerationFunction, postTimestepFunction, otherActionsFunction,
                ionStartMonitoringFunction,
                collisionModel);
        AppUtils::SignalHandler::setReceiver(verletIntegrator);
        verletIntegrator.run(timeSteps, dt);
    }
    else if (integratorMode==AppUtils::FULL_SUM_VERLET) {
        Integration::FullSumVerletIntegrator verletIntegrator(
                particles,
                singleStepAccelerationFunction, postTimestepFunction, otherActionsFunction,
                ionStartMonitoringFunction,
                collisionModel);
        AppUtils::SignalHandler::setReceiver(verletIntegrator);
        verletIntegrator.run(timeSteps, dt);
    }
    else if (integratorMode==AppUtils::PARALLEL_RUNGE_KUTTA4) {
        Integration::ParallelRK4Integrator rk4Integrator(
                particles,
                multiStepAccelerationFunction, multiStepSpaceChargeAccelerationFunction,
                postTimestepFunction, otherActionsFunction,
                ionStartMonitoringFunction,
                collisionModel);
        AppUtils::SignalHandler::setReceiver(rk4Integrator);
        rk4Integrator.run(timeSteps, dt);
    }
    else if (integratorMode==AppUtils::FULL_SUM_RUNGE_KUTTA4) {
        Integration::FullSumRK4Integrator fullSumRK4Integrator(
                particles,
                multiStepAccelerationFunction, multiStepSpaceChargeAccelerationFunction,
                postTimestepFunction, otherActionsFunction,
                ionStartMonitoringFunction,
                collisionModel);
        AppUtils::SignalHandler::setReceiver(fullSumRK4Integrator);
        fullSumRK4Integrator.run(timeSteps, dt);
    }
#ifdef WITH_FMM_3d
    else if (integratorMode==AppUtils::FMM3D_VERLET) {
        Integration::FMMVerletIntegrator<FMM3D::FMMSolver> integrator(
                particles,
                singleStepAccelerationFunction, postTimestepFunction,
                otherActionsFunction, ionStartMonitoringFunction,
                collisionModel);

        if (simConf->isParameter("FMM3D_precision")) {
            integrator.getFMMSolver()->setRequestedPrecision(simConf->doubleParameter("FMM3D_precision"));
        }

        AppUtils::SignalHandler::setReceiver(integrator);
        integrator.run(timeSteps, dt);
    }
#endif
#ifdef WITH_EXAFMMT
    else if (integratorMode==AppUtils::EXAFMM_VERLET) {
        Integration::FMMVerletIntegrator<ExaFMMt::FMMSolver> integrator(
                particles,
                singleStepAccelerationFunction, postTimestepFunction,
                otherActionsFunction, ionStartMonitoringFunction,
                collisionModel);

        if (simConf->isParameter("ExaFMM_order")) {
            integrator.getFMMSolver()->setExpansionOrder(simConf->intParameter("ExaFMM_order"));
        }

        AppUtils::SignalHandler::setReceiver(integrator);
        integrator.run(timeSteps, dt);
    }
#endif
}