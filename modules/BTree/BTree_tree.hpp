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

 BTree_tree.hpp

 Non parallelized version of a Barnes-Hut tree

 The tree consists of a root node and a tracking infrastructure which allows direct
 access to a particle and the leaf tree node which contains that particle by a particle
 index. This tracking is required for keeping the basic ability to synchronize particles in an
 external simulation with the particles in the tree.
 The tree is also able to be dynamically updated when particle positions change.

 ****************************/

#ifndef BTree_tree_hpp
#define BTree_tree_hpp

#include "BTree_node.hpp"
#include <memory>
#include <list>
#include <unordered_map>

//forward declare own classes:
namespace Core{
    class Vector;
}

namespace BTree {
    class Particle;

    class Tree {

    public:
        Tree(Core::Vector min, Core::Vector max);

        //simple getter:
        Node* getRoot();
        std::list<Particle*>* getParticleList();
        int getNumberOfParticles();

        //charge calculation methods:
        void computeChargeDistribution();
        Core::Vector computeEFieldFromTree(Particle &particle);

        //particle modification methods:
        void insertParticle(Particle &particle, int ext_index);
        void removeParticle(int ext_index);
        Particle* getParticle(int ext_index);
        void updateParticleLocation(int ext_index, Core::Vector newLocation);
        
        void printParticles();

    private:
        std::unique_ptr<Node> root_; ///< the root node of the tree
        std::unique_ptr<std::list<Particle*>> iVec_; ///< a linked particle list, stores the particles in a linear order
        std::unique_ptr<std::unordered_map<int, std::list<Particle *>::const_iterator>> iMap_;
        ///< a map between ion indices (keys used by an external simulation) and pointers into the internal particle list
    };
}

#endif /* BTree_tree_hpp */