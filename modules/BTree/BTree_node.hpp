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
 BTree_node.hpp

 Node in the non parallelized Barnes-Hut tree
 ****************************/

#ifndef BTree_node_hpp
#define BTree_node_hpp

#include "Core_vector.hpp"
#include "BTree_genericBaseNode.hpp"
#include <ostream>

namespace BTree {

    class Particle;

    class Node: public GenericBaseNode<Node> {

    public:

        // Constructors and assignment operators:
        Node(Core::Vector min, Core::Vector max, Node* parent);
        Node(const Node& that) = delete;
        Node& operator=(const Node& that) = delete;

        // Member methods:
        Core::Vector computeElectricFieldFromTree(Core::Particle &targetP);
        
        //mostly diagnostic methods:
        void testSpatialTreeIntegrity();
        void testNodeParticleIntegrity();
        bool isNodeInSubtree(const Node *nodeToFind, bool debug) const;
        bool isParticleInSubtree(const TreeParticle *particle, bool debug) const;
    };
}

#endif /* BTree_node_hpp */
