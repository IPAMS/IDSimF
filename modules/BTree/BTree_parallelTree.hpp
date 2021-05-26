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
        [[nodiscard]] int getNumberOfParticles() const;

        int init();
        std::vector<int> countNodesOnLevels();
        Core::Vector computeEFieldFromTree(Particle &particle);

        void insertParticle(Particle &particle, int ext_index);
        void removeParticle(int ext_index);
        Particle* getParticle(int ext_index) const;
        void updateParticleLocation(int extIndex, Core::Vector newLocation, int* numNodesChanged);
        int updateNodes(int ver);
        
        void printParticles() const;

    private:
        std::unique_ptr<ParallelNode> root_;
        std::unique_ptr<std::list<Particle*>> iVec_; ///< a linked particle list, stores the particles in a linear order
        std::unique_ptr<std::unordered_map<int, std::list<Particle *>::const_iterator>> iMap_; ///< a map between the ion indices (keys used by SIMION) and the pointers into the internal particle list

        std::vector<int> nodesOnLevels_;
        std::vector<int> nodeStartIndicesOnLevels_; ///< Vector of serialized indices of the first nodes on the individual tree levels
        int numberOfNodesTotal_ = 0;
        int nodesReserveSize_ = 0; ///< Number of nodes which should be reserved in the serialized field calculation
        std::vector<BTree::ParallelNode*>nodesSerialized_; ///< Serial vector of links to the tree nodes, used for parallelized access to the nodes
        int nTreeLevels_ = 0; //number of levels in the tree

        [[nodiscard]] int getTreeDepth_() const;
        void updateLevelStartIndices_();
        void serializeNodes_();
        void updateNodeChargeState_();
    };
    
}

#endif /* BTree_parallelTree_hpp */
