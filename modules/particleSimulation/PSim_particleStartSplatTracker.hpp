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
 BTree_particleStartSplatTracker.hpp

 Tracker system to manage / record particle start and stop (splat) times / positions

 ****************************/

#ifndef PSim_ionStartSplatTracker_hpp
#define PSim_ionStartSplatTracker_hpp

#include "Core_vector.hpp"
#include <map>


//forward declare own classes:
namespace BTree{
    class Particle;
}

namespace ParticleSimulation{

    struct pMapEntry {
        double startTime;            ///< The start time of a particle
        double splatTime;            ///< The splat time of a particle
        Core::Vector startLocation;  ///< Start location of a particle
        Core::Vector splatLocation;  ///< Splat Location of a Particle
    };

    class ParticleStartSplatTracker {

    public:

        void particleStart(BTree::Particle* particle, double time);
        void particleSplat(BTree::Particle* particle, double time);

        pMapEntry get(BTree::Particle* particle);

    private:

        std::map<BTree::Particle*, pMapEntry> pMap_;
    };
}

#endif //PSim_ionStartSplatTracker_hpp
