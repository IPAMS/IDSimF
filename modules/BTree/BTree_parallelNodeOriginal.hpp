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


#ifndef BTree_parallel_nodeOriginal_hpp
#define BTree_parallel_nodeOriginal_hpp

#include "Core_vector.hpp"
#include "BTree_abstractNode.hpp"
#include "BTree_particle.hpp"
#include <ostream>
#include <vector>

namespace BTree {

    class ParallelNodeOriginal: public AbstractNode {

    public:
        /**
         * Sub-Octants of a node (e.g. SWT = south west top, NEB = north east bottom)
        */
        enum Octant {SWT,NWT,SET,NET,SWB,NWB,SEB,NEB};

        // Conststructor:
        ParallelNodeOriginal(Core::Vector min, Core::Vector max, ParallelNodeOriginal* parent);
        
        // Copy Constructors:
        ParallelNodeOriginal(const ParallelNodeOriginal& that);
        ParallelNodeOriginal& operator=(const ParallelNodeOriginal& that);
        
        // Destructor:
        virtual ~ParallelNodeOriginal() override;

        // Simple accessors:
        int max_rek_tiefe(BTree::ParallelNodeOriginal* nod);
        void bestimme_Anfang(std::vector<int> &begin, std::vector<int> &zaehler, int n);
        void anz_Knoten_auf_der_Ebene(BTree::ParallelNodeOriginal* node, int ebene, std::vector<int> &zaehler, int n);
        void insertNodes(BTree::ParallelNodeOriginal* nod, std::vector<BTree::ParallelNodeOriginal*> &Nodes, int ebene, std::vector<int> &zahler, std::vector<int> &begin, int n);
        void update(std::vector<BTree::ParallelNodeOriginal*> &Nodes, std::vector<int> &begin, int n, int summe);
        void setDurchmesser();
        
        // Member methods:
        void initAsRoot() override;
        bool isRoot() const override;
        Octant getOctant(Core::Vector location);
        ParallelNodeOriginal* createOctNode(Octant oct);
        void insertParticle(Particle* particle) override;
        void removeMyselfFromTree() override;
        void updateSelf() override;
        void updateParents() override;

        void computeChargeDistributionRecursive() override;
        //Vector computeE(Particle& targetP,Vector efield);
        Core::Vector computeEFieldFromTree(Particle& targetP, std::vector<BTree::ParallelNodeOriginal*> &MyNodes);
        Core::Vector computeElectricForceFromTree(Particle &targetP);
        
        //mostly diagnostic methods:
        std::string toString() const override;
        void printTree(int level) const;
        void writeToStream(std::ostream& filestream, void (*writeFct)(std::ostream& filestream, const BTree::ParallelNodeOriginal* node)) const;

    private:
        //member fields:
        double durch_;

        ParallelNodeOriginal* octNodes_[8]; ///< Subnodes of this node
        ParallelNodeOriginal* parent_; ///< Parent node of this node

    };
}



#endif /* BTree_parallel_nodeOriginal_hpp */
