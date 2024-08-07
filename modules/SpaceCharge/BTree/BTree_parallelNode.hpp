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

 BTree_parallelNode.hpp

 Node in the parallel version of the Barnes-Hut tree

 @author Miguel Gruse
 @author Walter Wissdorf
 ****************************/


#ifndef BTree_parallel_node_hpp
#define BTree_parallel_node_hpp

#include "Core_vector.hpp"
#include "BTree_genericBaseNode.hpp"
#include "Core_particle.hpp"
#include <ostream>
#include <vector>

namespace BTree {

    class ParallelNode: public GenericBaseNode<ParallelNode> {

    friend class ParallelTree;

    public:

        // Constructors and assignment operators:
        ParallelNode(Core::Vector min, Core::Vector max, ParallelNode* parent);
        ParallelNode(const ParallelNode& that) = delete;
        ParallelNode& operator=(const ParallelNode& that) = delete;

        // Member methods:
        void serializeIntoVector(std::vector<BTree::ParallelNode*> &serializedNodes, std::size_t treeLevel,
                                 std::vector<std::size_t> &insertPositions);

        [[nodiscard]] std::size_t maximumRecursionDepth() const;
        void countNodesOnLevel(std::size_t level, std::vector<std::size_t> &numOfNodesOnLevels) const;

    private:

        double edgeLengthSquaredNormalized_ = 0.0;  ///< Squared edge length normalized by theta
        void updateNormalizedEdgeLength_();
    };
}



#endif /* BTree_parallel_node_hpp */
