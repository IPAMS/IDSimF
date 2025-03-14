/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2025 - Physical and Theoretical Chemistry /
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
 CollisionModel_AtomInteractionMap.hpp

 Simple map to map from pairs of Atoms to arbitrary values, mostly for storing pairwise atom interaction parameters

 ****************************/

#ifndef COLLISIONMODEL_ATOMINTERACTIONMAP_HPP
#define COLLISIONMODEL_ATOMINTERACTIONMAP_HPP

#include "CollisionModel_Atom.hpp"
#include <unordered_map>
#include <iostream>

namespace CollisionModel{
    template <typename valueType> class AtomInteractionMap {
    public:
        void insert(const Atom &atomA, const Atom &atomB, valueType value);
        valueType get(const Atom &atomA, const Atom &atomB);
        int size() const;
        //valueType operator[](const Atom &atomA, const Atom &atomB);

    private:
        std::unordered_map<std::size_t, std::unordered_map<std::size_t, double>> valueMap_;

        static void orderAtomPointers_(const Atom &atomA, const Atom &atomB, std::size_t& indexAtomL, std::size_t& indexAtomU);
    };

    template<typename valueType>
    void AtomInteractionMap<valueType>::orderAtomPointers_(const Atom& atomA, const Atom& atomB, std::size_t& indexAtomL, std::size_t& indexAtomU) {
        std::size_t indexA = atomA.getSpeciesIndex();
        std::size_t indexB = atomB.getSpeciesIndex();
        if (indexA < indexB) {
            indexAtomL = indexA;
            indexAtomU = indexB;
        }
        else {
            indexAtomL = indexB;
            indexAtomU = indexA;
        }
    }

    template<typename valueType>
    void CollisionModel::AtomInteractionMap<valueType>::insert(const Atom &atomA, const Atom &atomB, valueType value) {
        // we use the lower species index first "dimension" of the mapping
        std::size_t indexAtomL;
        std::size_t indexAtomU;
        orderAtomPointers_(atomA, atomB, indexAtomL, indexAtomU);

        if (auto valueL = valueMap_.find(indexAtomL); valueL != valueMap_.end()) {
            valueL->second.insert({indexAtomU, value}); //note: re assignments are silently ignored!!
        }
        else {
            valueMap_.emplace(indexAtomL, 0);
            valueMap_.at(indexAtomL).insert({indexAtomU, value});
        }
    }

    template<typename valueType>
    valueType CollisionModel::AtomInteractionMap<valueType>::get(const Atom &atomA, const Atom &atomB) {
        std::size_t indexAtomL;
        std::size_t indexAtomU;
        orderAtomPointers_(atomA, atomB, indexAtomL, indexAtomU);
        return valueMap_.at(indexAtomL).at(indexAtomU);
    }

    template<typename valueType>
    int CollisionModel::AtomInteractionMap<valueType>::size() const{
        return valueMap_.size();
    }
}

#endif //COLLISIONMODEL_ATOMINTERACTIONMAP_HPP
