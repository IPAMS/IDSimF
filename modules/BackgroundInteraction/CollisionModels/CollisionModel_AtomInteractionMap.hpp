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

namespace CollisionModel{
    template <typename valueType> class AtomInteractionMap {
    public:
        void insert(const Atom &atomA, const Atom &atomB, valueType value);
        valueType get(const Atom &atomA, const Atom &atomB);
        //valueType operator[](const Atom &atomA, const Atom &atomB);

    private:
        std::unordered_map<const Atom*, std::unordered_map<const Atom*, double>> valueMap_;

        static void orderAtomPointers_(const Atom &atomA, const Atom &atomB, const Atom* &pAtomL, const Atom* &pAtomU);
    };

    template<typename valueType>
    void AtomInteractionMap<valueType>::orderAtomPointers_(const Atom& atomA, const Atom& atomB, const Atom* &pAtomL, const Atom* &pAtomU) {
        if (&atomA < &atomB) {
            pAtomL = &atomA;
            pAtomU = &atomB;
        }
        else {
            pAtomL = &atomB;
            pAtomU = &atomA;
        }
    }

    template<typename valueType>
    void CollisionModel::AtomInteractionMap<valueType>::insert(const Atom &atomA, const Atom &atomB, valueType value) {
        // we use the atom pointer with the lower value as first "dimension" of the
        const Atom* pAtomL;
        const Atom* pAtomU;
        orderAtomPointers_(atomA, atomB, pAtomL, pAtomU);

        if (auto valueL = valueMap_.find(pAtomL); valueL != valueMap_.end()) {
            valueL->second.insert({pAtomU, value}); //note: re assignments are silently ignored!!
        }
        else {
            valueMap_.emplace(pAtomL, 0);
            valueMap_.at(pAtomL).insert({pAtomU, value});
        }
    }

    template<typename valueType>
    valueType CollisionModel::AtomInteractionMap<valueType>::get(const Atom &atomA, const Atom &atomB) {
        const Atom* pAtomL;
        const Atom* pAtomU;
        orderAtomPointers_(atomA, atomB, pAtomL, pAtomU);
        return valueMap_.at(pAtomL).at(pAtomU);
    }
}

#endif //COLLISIONMODEL_ATOMINTERACTIONMAP_HPP
