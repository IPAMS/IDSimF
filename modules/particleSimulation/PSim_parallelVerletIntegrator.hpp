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
 BTree_parallelVerletIntegrator.hpp
 
 Parallel version of a Barnes Hut tree based ion trajectory integration
 (simple verlet trajectory integrator) including space charge

 ****************************/

//FIXME: General: Establish dedicated types for time and use size_t

#ifndef BTree_parallelVerletIntegrator_hpp
#define BTree_parallelVerletIntegrator_hpp

#include <vector>
#include "BTree_particle.hpp"
#include "Core_vector.hpp"
#include "BTree_parallelTree.hpp"
#include "PSim_abstractTimeIntegrator.hpp"
#include "PSim_simpleVTKwriter.hpp"
#include "PSim_trajectoryExplorerJSONwriter.hpp"
#include "PSim_averageChargePositionWriter.hpp"
#include "CollisionModel_AbstractCollisionModel.hpp"

namespace ParticleSimulation{

    //std::function<Core::Vector(Core::Particle* particle, int particleIndex, Core::Tree& tree, double time, int timestep)> accelerationFctType;

    class ParallelVerletIntegrator: public AbstractTimeIntegrator {

        public:
            typedef std::function
                    <Core::Vector (BTree::Particle* particle,
                                   int particleIndex,
                                   BTree::ParallelTree& tree,
                                   double time,
                                   int timestep)>
                    accelerationFctType;

            /**
             * type definition for functions exporting data in every timestep
             */
            typedef std::function
                    <void (std::vector<BTree::Particle*>& particles,
                           BTree::ParallelTree& tree,
                           double time,
                           int timestep,
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
                           BTree::ParallelTree& tree,
                           double time,
                           int timestep)>
                    otherActionsFctType;

            ParallelVerletIntegrator(
                    std::vector<BTree::Particle*> particles,
                    accelerationFctType accelerationFunction,
                    timestepWriteFctType timestepWriteFunction = nullptr,
                    otherActionsFctType otherActionsFunction = nullptr,
                    AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction = nullptr,
                    CollisionModel::AbstractCollisionModel* collisionModel = nullptr
            );

            ParallelVerletIntegrator(
                    accelerationFctType accelerationFunction,
                    timestepWriteFctType timestepWriteFunction = nullptr,
                    otherActionsFctType otherActionsFunction = nullptr,
                    AbstractTimeIntegrator::particleStartMonitoringFctType ionStartMonitoringFunction = nullptr,
                    CollisionModel::AbstractCollisionModel* collisionModel = nullptr
            );

            void addParticle(BTree::Particle* particle) override;
            void run(int nTimesteps, double dt) override;
            void runSingleStep(double dt) override;
            void finalizeSimulation() override;

    private:

        size_t numberOfNodes_; ///< number of nodes in the parallel BTree attached to this integrator

        CollisionModel::AbstractCollisionModel* collisionModel_; ///< the gas collision model to perform while integrating

        accelerationFctType accelerationFunction_ = nullptr;   ///< function to calculate particle acceleration
        timestepWriteFctType timestepWriteFunction_ = nullptr; ///< function to export / write time step results
        otherActionsFctType otherActionsFunction_ = nullptr;   ///< function for arbitrary other actions in the simulation

        //internal variables for actual calculations:
        Core::Vector loc_min_ = Core::Vector(-1000,-1000,-1000); ///< currently hard coded lower corner of the sim. domain
        Core::Vector loc_max_ = Core::Vector( 1000, 1000, 1000); ///< currently hard coded upper corner of the sim. domain
        BTree::ParallelTree tree_ = BTree::ParallelTree(loc_min_,loc_max_);
        ///< The parallel BTree (primarily for space charge calculation)

        std::vector<Core::Vector>  newPos_;  ///< new position (after time step) for particles
        std::vector<Core::Vector>  a_t_;     ///< last time step acceleration for particles
        std::vector<Core::Vector>  a_tdt_;   ///< new acceleration for particles

        void bearParticles_(double time);
        void initInternalState_();
    };
}


#endif /* BTree_parallelVerletIntegrator_hpp */
