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
 BTree_abstractTimeIntegrator.hpp

 Simple base class for time trajectory integrators

 ****************************/

#ifndef BTree_abstractTimeIntegrator_hpp
#define BTree_abstractTimeIntegrator_hpp

#include "BTree_particle.hpp"
#include "Core_vector.hpp"
#include <cstdio>
#include <vector>
#include <functional>

namespace ParticleSimulation{

    /**
     * Abstract base class for trajectory time integrators
     */
    class AbstractTimeIntegrator {

    public:
        enum RunState {RUNNING, STOPPED, IN_TERMINATION};
        typedef std::pair<double, BTree::Particle*> pTobPair_t;

        /**
          * Type definition for a function to watch ion start events
          */
        typedef std::function
                <void (BTree::Particle* particle,
                       double time)>
                particleStartMonitoringFctType;

        explicit AbstractTimeIntegrator(particleStartMonitoringFctType ionStartMonitorFct = nullptr);
        explicit AbstractTimeIntegrator(const std::vector<BTree::Particle*>& particles,
                               particleStartMonitoringFctType ionStartMonitorFct = nullptr);

        virtual ~AbstractTimeIntegrator() = default;

        virtual void addParticle(BTree::Particle* particle) = 0;
        virtual void run(int nTimesteps, double dt) = 0;
        virtual void runSingleStep(double dt) = 0;
        virtual void finalizeSimulation() = 0;

        void setTerminationState();
        [[nodiscard]] RunState runState() const;
        [[nodiscard]] double time() const;
        [[nodiscard]] int timeStep() const;

    protected:
        RunState runState_ = STOPPED; ///< the current state the integrator is in
        double time_ = 0.0; ///< the current time in the simulation
        int timestep_ = 0; ///< the current time step
        std::vector<BTree::Particle*> particles_; ///< links to the simulated particles
        std::size_t nParticles_ = 0; ///< number of particles
        std::vector<pTobPair_t> particleTOBs_; ///< Time of births of the individual particles
        std::size_t particlesBornIdx_ = 0; ///< index in particleTOBs_ indicating the particles already born
        particleStartMonitoringFctType particleStartMonitorFct_ = nullptr; ///< Monitoring function for
        bool bearParticles_(double time);
    };
}
#endif /* BTree_abstractTimeIntegrator_hpp */