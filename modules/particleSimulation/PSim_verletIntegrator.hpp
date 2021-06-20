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

#include "PSim_abstractTimeIntegrator.hpp"
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

namespace ParticleSimulation{

    //std::function<Core::Vector(Core::Particle* particle, int particleIndex, Core::Tree& tree, double time, int timestep)> accelerationFctType;

    /**
     * Core based ion trajectory integration with a simple verlet trajectory integrator scheme including space charge.
     *
     * The acceleration calculation, additional ("other") actions performed in every time step and what is exported
     * in every time step is passed to this trajectory integrator externally by functions. Thus, the integration scheme
     * can be applied to arbitrary simulation problems.
     */
    class VerletIntegrator: public AbstractTimeIntegrator{

    public:

        /**
        * type definition for acceleration calculation functions
        */
        typedef std::function
                <Core::Vector (BTree::Particle* particle,
                                int particleIndex,
                                BTree::Tree& tree,
                                double time,
                                unsigned int timestep)>
                accelerationFctType;

        /**
         * type definition for functions exporting data in every timestep
         */
        typedef std::function
                <void (std::vector<BTree::Particle*>& particles,
                       BTree::Tree& tree,
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
                       int particleIndex,
                       BTree::Tree& tree,
                       double time,
                       unsigned int timestep)>
                otherActionsFctType;

            VerletIntegrator(
                    std::vector<BTree::Particle*> particles,
                    accelerationFctType accelerationFunction,
                    timestepWriteFctType timestepWriteFunction = nullptr,
                    otherActionsFctType otherActionsFunction = nullptr,
                    AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction = nullptr,
                    CollisionModel::AbstractCollisionModel* collisionModel = nullptr
            );

            VerletIntegrator(
                    accelerationFctType accelerationFunction,
                    timestepWriteFctType timestepWriteFunction = nullptr,
                    otherActionsFctType otherActionsFunction = nullptr,
                    AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction = nullptr,
                    CollisionModel::AbstractCollisionModel* collisionModel = nullptr
            );

            void addParticle(BTree::Particle* particle) override;
            void run(unsigned int nTimesteps, double dt) override;
            void runSingleStep(double dt) override;
            void finalizeSimulation() override;

    private:
        CollisionModel::AbstractCollisionModel* collisionModel_ = nullptr; ///< a gas collision model active in the simulation

        accelerationFctType accelerationFunction_ = nullptr; ///< function to calculate particle acceleration
        timestepWriteFctType timestepWriteFunction_ = nullptr; ///< function to define what is exported in every time step
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