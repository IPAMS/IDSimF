/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2021 - Physical and Theoretical Chemistry /
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
 appUtils_commandlineParsing.hpp

 Parsing of commandline arguments for app runs

 ****************************/

#ifndef IDSIMF_APPUTILS_COMMANDLINEPARSER_HPP
#define IDSIMF_APPUTILS_COMMANDLINEPARSER_HPP

#include "appUtils_simulationConfiguration.hpp"

namespace AppUtils{
    class CommandlineParser {
    public:
        CommandlineParser(int argc, const char * argv[], std::string appName, std::string appDescription,  bool multithreaded=false);
        AppUtils::SimulationConfiguration simulationConfiguration();
        std::string projectName();
        std::string confFileName();
        unsigned int numberOfThreads();

    private:
        AppUtils::SimulationConfiguration simulationConfiguration_;
        std::string projectName_;
        std::string confFileName_;
        unsigned int numberOfThreads_= 0;
    };
}




#endif //IDSIMF_APPUTILS_COMMANDLINEPARSER_HPP
