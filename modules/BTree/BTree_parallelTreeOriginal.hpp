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

#ifndef BTree_parallelTreeOriginal_hpp
#define BTree_parallelTreeOriginal_hpp

#include <stdio.h>
#include <list>
#include <unordered_map>
#include "BTree_parallelNodeOriginal.hpp"

namespace BTree {

    //TODO: Adding particles with exactly the same position leads to crash in verlet integrator, (throw exception here??)

    class ParallelTreeOriginal {
    private:
        std::unique_ptr<ParallelNodeOriginal> root_;
        std::unique_ptr<std::list<Particle*>> iVec_; ///< a linked particle list, stores the particles in a linear order
        std::unique_ptr<std::unordered_map<int, std::list<Particle *>::const_iterator>> iMap_; ///< a map between the ion indices (keys used by SIMION) and the pointers into the internal particle list
    public:
        // Conststructor:
        ParallelTreeOriginal(Core::Vector min, Core::Vector max);

        // Public member methods:
        ParallelNodeOriginal* getRoot();
        void computeChargeDistribution();
        Core::Vector computeEFieldFromTree(Particle &particle, std::vector<BTree::ParallelNodeOriginal*> &MyNodes);

        std::list<Particle*>* getParticleList();
        int getNumberOfParticles();
        
        void bestimme_max_rek_tiefe(int* n);
        void bestimme_anzahl_Knoten_auf_ebene(int ebene, std::vector<int> &zaehler, int n);
        void bestimme_Feld_begin(std::vector<int> &begin, std::vector<int> &zaehler, int n);
        void einfuegen(std::vector<BTree::ParallelNodeOriginal*> &Nodes, int ebene, std::vector<int> &zaehler, std::vector<int> &begin, int n);
        void Knotenup(std::vector<BTree::ParallelNodeOriginal*> &Nodes, std::vector<int> &begin, int n, int summe);
        void updateS(int ext_index, Core::Vector newLocation);
        void insertParticle(Particle &particle, int ext_index);
        void removeParticle(int ext_index);
        Particle* getParticle(int ext_index);
        void updateParticleLocationnew(int ext_index,Core::Vector newLocation,int* ver);//todo
        
        void printParticles();
    };
    
}

#endif /* BTree_parallelTreeOriginal_hpp */
