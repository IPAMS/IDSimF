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
 RS_ConfigFileParser.hpp

 This implements a configuration file parser / simulation configuration generator for
 chemical (reaction simulation) simulations

 ****************************/

#ifndef RS_ConfigFileParser_hpp
#define RS_ConfigFileParser_hpp

#include "RS_Substance.hpp"
#include "RS_SimulationConfiguration.hpp"
#include "RS_AbstractReaction.hpp"
#include "RS_StaticReaction.hpp"
#include "RS_StaticThermalizingReaction.hpp"
#include "RS_SimpleCollisionStepReaction.hpp"
#include "RS_VantHoffReaction.hpp"
#include "RS_FieldDependentVantHoffReaction.hpp"
#include <iostream>
#include <memory>
#include <vector>
#include <utility>
#include <fstream>
#include <sstream>
#include <regex>
#include <algorithm>


namespace RS {
    // individual exception for problems in ion cloud files
    class ConfigurationFileException : public std::exception
    {
    public:
        explicit ConfigurationFileException (std::string s) {this->message = s;}
        ~ConfigurationFileException() override = default;
        const char * what() const throw() override {return message.c_str();}

    private:
        std::string message;
    };

    class ConfigFileParser {

    public:
        std::pair<std::vector<std::string>,std::vector<std::string>>
            splitString(std::string str, const std::string& patterntxt) const;
        std::unique_ptr<SimulationConfiguration> parseText(const std::string& confStr) const;
        std::unique_ptr<SimulationConfiguration> parseFile(const std::string& filename) const;

        std::unique_ptr<SimulationConfiguration> getTestConfigWaterClusters() const;
        std::unique_ptr<SimulationConfiguration> getTestConfigSimple() const;

    private:
        bool parseReactions(SimulationConfiguration* simConf, const std::string& input) const;
        std::pair<int,std::string> parseReactionPartnerString(const std::string& input) const;
        bool parseSubstances(SimulationConfiguration* simConf, const std::string& input)const;
    };
}


#endif /*RS_ConfigFileParser_hpp*/
