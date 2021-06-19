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

 BTree_parallelTree.hpp

    Parallel version of the Barnes-Hut tree
 ****************************/

#ifndef BTree_parallelTree_hpp
#define BTree_parallelTree_hpp

#include <stdio.h>
#include <list>
#include <unordered_map>
#include "BTree_parallelNode.hpp"

namespace BTree {

    class ParallelTree {
    public:
        // Constructor:
        ParallelTree(Core::Vector min, Core::Vector max);

        // Public member methods:
        [[nodiscard]] ParallelNode* getRoot() const;
        std::list<Particle*>* getParticleList() const;
        [[nodiscard]] std::size_t getNumberOfParticles() const;

        int init();
        std::vector<std::size_t> countNodesOnLevels();
        Core::Vector computeEFieldFromTree(Particle &particle);

        void insertParticle(Particle &particle, std::size_t ext_index);
        void removeParticle(std::size_t ext_index);
        Particle* getParticle(std::size_t ext_index) const;
        void updateParticleLocation(std::size_t extIndex, Core::Vector newLocation, int* numNodesChanged);
        int updateNodes(int ver);
        
        void printParticles() const;

    private:
        std::unique_ptr<ParallelNode> root_;
        std::unique_ptr<std::list<Particle*>> iVec_; ///< a linked particle list, stores the particles in a linear order
        std::unique_ptr<std::unordered_map<std::size_t, std::list<Particle *>::const_iterator>> iMap_; ///< a map between the ion indices (keys used by SIMION) and the pointers into the internal particle list

        std::vector<std::size_t> nodesOnLevels_;
        std::vector<std::size_t> nodeStartIndicesOnLevels_; ///< Vector of serialized indices of the first nodes on the individual tree levels
        std::size_t numberOfNodesTotal_ = 0;
        std::size_t nodesReserveSize_ = 0; ///< Number of nodes which should be reserved in the serialized field calculation
        std::vector<BTree::ParallelNode*>nodesSerialized_; ///< Serial vector of links to the tree nodes, used for parallelized access to the nodes
        std::size_t nTreeLevels_ = 0; //number of levels in the tree

        [[nodiscard]] std::size_t getTreeDepth_() const;
        void updateLevelStartIndices_();
        void serializeNodes_();
        void updateNodeChargeState_();
    };
    
}

#endif /* BTree_parallelTree_hpp */
