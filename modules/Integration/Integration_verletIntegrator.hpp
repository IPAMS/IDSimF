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
 BTree_verletIntegrator.hpp

 Barnes Hut Tree based ion trajectory integration with a simple verlet trajectory integrator scheme
 including space charge

 ****************************/

#ifndef BTree_verletIntegrator_hpp
#define BTree_verletIntegrator_hpp

#include "Integration_abstractTimeIntegrator.hpp"
#include "Core_vector.hpp"
#include "BTree_tree.hpp"

#include <cstdio>
#include <vector>
#include <functional>
//forward declare own classes:
namespace BTree{
    class Particle;
}

namespace CollisionModel{
    class AbstractCollisionModel;
}

namespace Integration{
    /**
     * Ion trajectory integrator with a simple velocity verlet trajectory integrator scheme including space charge.
     * Space charge calculation with a serial Barnes Hut algorithm.
     *
     * The acceleration calculation and additional actions performed is passed to this trajectory integrator externally
     * by functions. Thus, the integration scheme can be applied to arbitrary simulation problems.
     */
    class VerletIntegrator: public AbstractTimeIntegrator{

    public:

        VerletIntegrator(
                std::vector<Core::Particle*> particles,
                accelerationFctSingleStepType accelerationFunction,
                postTimestepFctType postTimestepFunction = nullptr,
                otherActionsFctType otherActionsFunction = nullptr,
                AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction = nullptr,
                CollisionModel::AbstractCollisionModel* collisionModel = nullptr
        );

        VerletIntegrator(
                accelerationFctSingleStepType accelerationFunction,
                postTimestepFctType postTimestepFunction = nullptr,
                otherActionsFctType otherActionsFunction = nullptr,
                AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction = nullptr,
                CollisionModel::AbstractCollisionModel* collisionModel = nullptr
        );

        void addParticle(Core::Particle* particle) override;
        void run(unsigned int nTimesteps, double dt) override;
        void runSingleStep(double dt) override;
        void finalizeSimulation() override;

    private:
        CollisionModel::AbstractCollisionModel* collisionModel_ = nullptr; ///< a gas collision model active in the simulation

        accelerationFctSingleStepType accelerationFunction_ = nullptr; ///< function to calculate particle acceleration
        postTimestepFctType postTimestepFunction_ = nullptr; ///< function to define what is exported in every time step
        otherActionsFctType otherActionsFunction_ = nullptr; ///< function with other actions done in every time step


        //internal variables for actual calculations:
        Core::Vector loc_min_ = Core::Vector(-1000,-1000,-1000);
        Core::Vector loc_max_ = Core::Vector( 1000, 1000, 1000);
        BTree::Tree tree_ = BTree::Tree(loc_min_,loc_max_); ///< tree used for particles and space charge calculation

        //new position, and intermediate results for every particle:
        std::vector<Core::Vector>  newPos_;
        std::vector<Core::Vector>  a_t_;
        std::vector<Core::Vector>  a_tdt_;

        //std::vector<int> ion_keys_; ///< vector of external labels attached to the individual ions

        void bearParticles_(double time);
    };
}
#endif /* BTree_verletIntegrator_hpp */