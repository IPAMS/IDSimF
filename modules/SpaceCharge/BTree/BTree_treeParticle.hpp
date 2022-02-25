/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2022 - Physical and Theoretical Chemistry /
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

 BTree_treeParticle.hpp

 Simulated particle in a BTree (mostly a wrapper for Core::Particle or RS::ReactiveParticle)

 ****************************/

#ifndef IDSIMF_BTREE_TREEPARTICLE_HPP
#define IDSIMF_BTREE_TREEPARTICLE_HPP

#include "Core_particle.hpp"

//Forward declare own classes
namespace BTree {
    class AbstractNode;
}

namespace BTree{
    /**
     * Defines a wrapper for a simulated particle, which can be managed by a Barnes Hut Tree (BTree)
     */
    class TreeParticle {

    public:
        TreeParticle(Core::Particle* baseParticle);

        //[[nodiscard]] Core::Particle* get() const;

        /** Gets the wrapped particle*/
        /*[[nodiscard]] inline Core::Particle* get() const {
            return particle_;
        }*/
        void setHostNode(BTree::AbstractNode* newHostNode);
        [[nodiscard]] BTree::AbstractNode* getHostNode() const;

        Core::Particle* wrappedParticle = nullptr;

    private:
        //Core::Particle* particle_ = nullptr; ////< The particle which is wrapped by this TreeParticle
        BTree::AbstractNode* hostNode_ = nullptr; ///< A link to a tree node to which the particle is belonging currently
    };
}

#endif //IDSIMF_BTREE_TREEPARTICLE_HPP
