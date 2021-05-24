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
 RS_Simulation.hpp

 This implements the RS simulation itself, which manages existing substances / elapsed time etc.

 ****************************/

#ifndef RS_Simulation_hpp
#define RS_Simulation_hpp

#include "RS_Substance.hpp"
#include "RS_AbstractReaction.hpp"
#include "RS_ReactiveParticle.hpp"
#include "RS_ConfigFileParser.hpp"
#include "spdlog/spdlog.h"
#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <list>
#include <utility>

namespace RS{class Simulation;}
std::ostream& operator<<(std::ostream& os, const RS::Simulation& sim);


namespace RS {
    class Simulation {

    using reactionMap = std::map<Substance* const, std::vector<AbstractReaction*>>;

    public:
        explicit Simulation(const std::string& configFileName);
        explicit Simulation(std::unique_ptr<RS::SimulationConfiguration> simConf);

        SimulationConfiguration* simulationConfiguration();

        bool addParticle(RS::ReactiveParticle* particle,int index);
        void removeParticle(int index);
        ReactiveParticle& getParticle(int index);
        std::map<Substance* const,int> discreteConcentrations();
        long totalReactionEvents();
        long illEvents();
        long reactionEvents(AbstractReaction* reaction) const;

        void performTimestep(ReactionConditions& conditions, double dt);
        void doReaction(RS::AbstractReaction* reaction, RS::ReactiveParticle* particle, RS::Substance* product);
        bool react(int index, ReactionConditions& conditions, double dt);
        bool collisionReact(int index, RS::Substance* reactionPartnerSpecies, CollisionConditions& conditions);
        void advanceTimestep(double dt);

        int timestep();
        double simulationTime();

        void printConcentrations();
        void printReactionStatistics();
        void logConcentrations(std::shared_ptr<spdlog::logger>& logger);
        void logReactionStatistics(std::shared_ptr<spdlog::logger>& logger);
        friend std::ostream& ::operator<<(std::ostream& os, const RS::Simulation& sim);


    private:

        using pMap = std::unordered_map<int, RS::ReactiveParticle*>;
        using pPair= pMap::value_type;
        pMap particleMap_;

        void initFromSimulationConfig_(std::unique_ptr<RS::SimulationConfiguration> simConf);
        reactionMap indReactDeepCopy_();
        bool react_(int index, ReactionConditions& conditions, double dt, reactionMap &reacInd);

        std::string concentrationString_();
        std::string reactionStatisticsString_();

        //implement private members / data structures / methods
        long totalReactionEvents_ = 0; ///< the total number of reaction events in the simulation
        long illEvents_= 0;  ///< the number of illegal / ill events (events with reacion probabilties > 1)
        double sumTime_ = 0.0; ///< the cumulative time (sum of all time steps)
        int nTimesteps_ = 0; ///< the total number of timesteps in the simulation
        std::unique_ptr<RS::SimulationConfiguration> simConf_;
        std::vector<Substance*> substances_; ///< the vector of all chemical substances in the simulation
        std::vector<AbstractReaction*> reactions_; ///< the vector of all chemical reactions in the simulation
        reactionMap reacInd_; ///< substance specific independent reactions (only one discrete educt)
        reactionMap reacDep_; ///< substance specific dependent reactions (more than one discrete educt)
        std::map<AbstractReaction* const,long> reactionEvents_; ///< a map to count the number of individual reaction events
        std::map<std::pair<Substance* const, Substance* const>, std::vector<AbstractReaction*>> reacCollision_; ///< map of collision based reactions
        std::map<Substance* const,std::vector<double>> staticProbabilities_; ///< substance specific static reaction probabilities
        std::map<Substance* const,int> discreteConcentrations_; ///< substance specific discrete particle concentrations

        //std::list<ReactiveParticle> particles_;
    };
}

#endif /* RS_Simulation_hpp */
