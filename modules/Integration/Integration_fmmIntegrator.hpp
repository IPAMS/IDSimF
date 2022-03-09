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

#ifndef IDSIMF_INTEGRATION_FMMINTEGRATOR_HPP
#define IDSIMF_INTEGRATION_FMMINTEGRATOR_HPP

#include "Core_particle.hpp"
#include "Core_vector.hpp"
#include "Integration_abstractTimeIntegrator.hpp"
#include "CollisionModel_AbstractCollisionModel.hpp"
#include <vector>

namespace Integration{

    //std::function<Core::Vector(Core::Particle* particle, int particleIndex, Core::Tree& tree, double time, int timestep)> accelerationFctType;

    template <class FMMSolverT>
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
        FMMSolverT solver_;

        CollisionModel::AbstractCollisionModel* collisionModel_ = nullptr; ///< the gas collision model to perform while integrating

        accelerationFctType accelerationFunction_ = nullptr;   ///< function to calculate particle acceleration
        timestepWriteFctType timestepWriteFunction_ = nullptr; ///< function to export / write time step results
        otherActionsFctType otherActionsFunction_ = nullptr;   ///< function for arbitrary other actions in the simulation

        std::vector<Core::Vector>  a_t_;     ///< last time step acceleration for particles
        std::vector<Core::Vector>  a_tdt_;   ///< new acceleration for particles

    };


    template <class FMMSolverT>
    FMMVerletIntegrator<FMMSolverT>::FMMVerletIntegrator(
            const std::vector<Core::Particle*>& particles,
            accelerationFctType accelerationFunction,
            timestepWriteFctType timestepWriteFunction,
            otherActionsFctType otherActionsFunction,
            AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction,
            CollisionModel::AbstractCollisionModel* collisionModel):
        AbstractTimeIntegrator(particles, ionStartMonitoringFunction),
        collisionModel_(collisionModel),
        accelerationFunction_(std::move(accelerationFunction)),
        timestepWriteFunction_(std::move(timestepWriteFunction)),
        otherActionsFunction_(std::move(otherActionsFunction))
    {}

    template <class FMMSolverT>
    FMMVerletIntegrator<FMMSolverT>::FMMVerletIntegrator(
            accelerationFctType accelerationFunction,
            timestepWriteFctType timestepWriteFunction,
            otherActionsFctType otherActionsFunction,
            AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction,
            CollisionModel::AbstractCollisionModel* collisionModel) :
        AbstractTimeIntegrator(ionStartMonitoringFunction),
        collisionModel_(collisionModel),
        accelerationFunction_(std::move(accelerationFunction)),
        timestepWriteFunction_(std::move(timestepWriteFunction)),
        otherActionsFunction_(std::move(otherActionsFunction))
    {}

    template <class FMMSolverT>
    void FMMVerletIntegrator<FMMSolverT>::addParticle(Core::Particle* particle) {
        particles_.push_back(particle);
        a_t_.emplace_back(Core::Vector(0,0,0));
        a_tdt_.emplace_back(Core::Vector(0,0,0));

        solver_.insertParticle(*particle, nParticles_);
        ++nParticles_;
    }

    template <class FMMSolverT>
    void FMMVerletIntegrator<FMMSolverT>::run(unsigned int nTimesteps, double dt) {
        // run init:
        this->runState_ = RUNNING;
        bearParticles_(0.0);

        if (timestepWriteFunction_ !=nullptr) {
            timestepWriteFunction_(particles_, time_, timestep_, false);
        }

        // run:
        for (unsigned int step=0; step< nTimesteps; step++){
            runSingleStep(dt);
            if (this->runState_ == IN_TERMINATION){
                break;
            }
        }
        this->finalizeSimulation();
        this->runState_ = STOPPED;
    }

    template <class FMMSolverT>
    void FMMVerletIntegrator<FMMSolverT>::runSingleStep(double dt) {
        bearParticles_(time_);

        // first: Calculate charge distribution with particle positions and charges at the beginning of the time step
        solver_.computeChargeDistribution();

        // then: perform particle update / velocity verlet time integration
        #pragma omp parallel \
                default(none) shared(a_tdt_, a_t_, dt, particles_)
        {

            #pragma omp for
            for (std::size_t i=0; i<nParticles_; i++){

                if (particles_[i]->isActive()){

                    if (collisionModel_ != nullptr) {
                        collisionModel_->updateModelParameters(*(particles_[i]));
                    }

                    Core::Vector newPos = particles_[i]->getLocation() + particles_[i]->getVelocity() * dt + a_t_[i]*(1.0/2.0*dt*dt);
                    a_tdt_[i] = accelerationFunction_(particles_[i], i, solver_, time_, timestep_);
                    //acceleration changes due to background interaction:

                    if (collisionModel_ != nullptr) {
                        collisionModel_->modifyAcceleration(a_tdt_[i], *(particles_[i]), dt);
                    }

                    particles_[i]->setVelocity( particles_[i]->getVelocity() + ((a_t_[i]+ a_tdt_[i])*1.0/2.0 *dt) );
                    a_t_[i] = a_tdt_[i];

                    //velocity changes due to background interaction:
                    if (collisionModel_ != nullptr) {
                        collisionModel_->modifyVelocity(*(particles_[i]),dt);
                        collisionModel_->modifyPosition(newPos, *(particles_[i]), dt);
                    }

                    if (otherActionsFunction_ != nullptr) {
                        otherActionsFunction_(newPos, particles_[i], i, time_, timestep_);
                    }

                    particles_[i]->setLocation(newPos);
                }
            }
        }

        // Update time and timestep, if existting: Call write timestep function
        time_ = time_ + dt;
        timestep_++;
        if (timestepWriteFunction_ != nullptr) {
            timestepWriteFunction_(particles_, time_, timestep_, false);
        }
    }

    template <class FMMSolverT>
    void FMMVerletIntegrator<FMMSolverT>::finalizeSimulation() {
        if (timestepWriteFunction_ != nullptr){
            timestepWriteFunction_(particles_, time_, timestep_, true);
        }
    }
}


#endif //IDSIMF_INTEGRATION_FMMINTEGRATOR_HPP
