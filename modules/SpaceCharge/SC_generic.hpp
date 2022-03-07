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
 SC-generic.hpp

 Generic interfaces for space charge calculation
 ****************************/


#ifndef IDSIMF_SC_GENERIC_HPP
#define IDSIMF_SC_GENERIC_HPP

#include "Core_vector.hpp"
#include "Core_particle.hpp"
#include <list>
#include <vector>
#include <unordered_map>

namespace Core{
    class Vector;
    class Particle;
}

namespace SpaceCharge{

    class FieldCalculator {

    public:
        virtual ~FieldCalculator() = default;
        virtual Core::Vector getEFieldFromSpaceCharge(Core::Particle &particle) = 0;
    };

    constexpr double NEGATIVE_ELECTRIC_CONSTANT = -1.0* Core::ELECTRIC_CONSTANT;

    struct particleListEntry{
        Core::Particle* particle;
        Core::Vector gradient;
        double potential;
    };

    using particlePtrList = std::list<particleListEntry>;
    /**
     * Abstract base class for Barnes-Hut Tree nodes
     */
    class GenericSpaceChargeSolver : public SpaceCharge::FieldCalculator {

    public:
        GenericSpaceChargeSolver();

        void insertParticle(Core::Particle &particle, std::size_t ext_index);
        void removeParticle(std::size_t ext_index);
        [[nodiscard]] std::size_t getNumberOfParticles() const;

        [[nodiscard]] virtual Core::Vector getEFieldFromSpaceCharge(Core::Particle& particle) override;
        virtual void computeChargeDistribution() = 0;

    protected:
        std::unique_ptr<particlePtrList> iVec_; ///< a linked particle list, stores the particles in a linear order
        std::unique_ptr<std::unordered_map<std::size_t, particlePtrList::const_iterator>> iMap_; ///< a map between the ion indices (keys used by SIMION) and the pointers into the internal particle list
        std::unique_ptr<std::unordered_map<Core::Particle*, particlePtrList::const_iterator>> pMap_; ///< a map between the particle pointers and the pointers into the internal particle list
    };
}

#endif //IDSIMF_SC_GENERIC_HPP
