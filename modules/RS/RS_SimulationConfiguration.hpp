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
 RS_SimulationConfiguration.hpp

 Configuration for a RS simulation, which contains the set of configured chemical species
 and chemical reactions

 ****************************/

#ifndef RS_SimulationConfiguration_hpp
#define RS_SimulationConfiguration_hpp

#include <vector>
#include <memory>
#include "RS_Substance.hpp"
#include "RS_AbstractReaction.hpp"

namespace RS {

    // individual exception for problems in simulation configurations
    class SimulationConfigurationException : public std::exception
    {
    public:
        explicit SimulationConfigurationException (std::string s) {this->message = s;}
        ~SimulationConfigurationException() override = default;
        const char * what() const throw() override {return message.c_str();}

    private:
        std::string message;
    };

    class SimulationConfiguration {

    public:
        bool addSubstance(std::unique_ptr<Substance>& subst);
        [[nodiscard]] Substance* substance(std::size_t index) const;
        [[nodiscard]] Substance* substanceByName(std::string substanceName) const;
        [[nodiscard]] std::vector<Substance*> getAllSubstances() const;
        [[nodiscard]] std::vector<Substance*> getAllDiscreteSubstances() const;

        bool addReaction(std::unique_ptr<AbstractReaction>& reac);
        [[nodiscard]] AbstractReaction* reaction(std::size_t index) const;
        [[nodiscard]] std::vector<AbstractReaction*> getAllReactions() const;

        void updateConfiguration();

    private:
        std::vector <std::unique_ptr<Substance>> substances_; ///< a vector of substances (owned by this simulation configuration)
        std::vector <Substance*> discreteSubstances_; ///< Vector of substances modeled as discrete particles
        std::vector <std::unique_ptr<AbstractReaction>> reactions_; ///< a vector of abstract reactions
        std::map <std::string,Substance*> substancesNameMap_; ////< a map from the substance names to the substance pointers
    };
}

#endif //RS_SimulationConfiguration_hpp
