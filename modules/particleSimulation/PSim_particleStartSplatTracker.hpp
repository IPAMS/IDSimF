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

 Tracker system to record particle start and stop (splat) times / positions

 ****************************/

#ifndef PSim_ionStartSplatTracker_hpp
#define PSim_ionStartSplatTracker_hpp

#include "Core_vector.hpp"
#include <unordered_map>
#include <vector>


//forward declare own classes:
namespace BTree{
    class Particle;
}

namespace ParticleSimulation{



    /**
     * Tracker system to record particle start and stop (splat) times / positions. This separate tracking is
     * required, because simulations are allowed to destroy and generate particles at their will, for example
     * for simulations with a constant particle inflow. Thus an independent, global, tracking of start and
     * stop locations and times is required.
     */
    class ParticleStartSplatTracker {

    public:

        enum particleState {
            STARTED,
            SPLATTED
        };

        struct pMapEntry {
            int globalIndex;             ///< A particle index to identify the particle globally
            particleState state;         ///< The current state of the tracked particle
            double startTime=0;          ///< The start time of a particle
            double splatTime=0;          ///< The splat time of a particle
            Core::Vector startLocation;  ///< Start location of a particle
            Core::Vector splatLocation;  ///< Splat Location of a Particle
        };


        ParticleStartSplatTracker();

        void particleStart(BTree::Particle* particle, double time);
        void particleSplat(BTree::Particle* particle, double time);

        pMapEntry get(BTree::Particle* particle);
        void sortStartSplatData();

        std::vector<pMapEntry> getStartSplatData();

        std::vector<double> getStartTimes();
        std::vector<double> getSplatTimes();
        std::vector<Core::Vector> getStartLocations();
        std::vector<Core::Vector> getSplatLocations();

        //std::vector

    private:

        std::unordered_map<BTree::Particle*, pMapEntry> pMap_;
        std::vector<pMapEntry> sortedParticleData_;
        int pInsertIndex_ = 0;

    };
}

#endif //PSim_ionStartSplatTracker_hpp
