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
#include "appUtils_logging.hpp"
#include "CLI11.hpp"
#include <memory>

namespace AppUtils{

    /**
     * "Exception" / Message to indicate that the app should be terminated during commandline
     * parsing process.
     */
    class TerminatedWhileCommandlineParsing : public std::exception
            {
            public:
                TerminatedWhileCommandlineParsing (int returnCode, std::string message):
                    returnCode_(returnCode),
                    message_(message)
                    {}
                ~TerminatedWhileCommandlineParsing() override = default;
                const char * what() const throw() override {return message_.c_str();}
                int returnCode() const {return returnCode_;}

            private:
                int returnCode_;
                std::string message_;
            };

    /**
     * Parser for the commandline arguments of calls of simulation applications
     */
    class CommandlineParser {
    public:
        CommandlineParser(int argc, const char * argv[], std::string appName, std::string appDescription,  bool multithreaded=false);
        simConf_ptr simulationConfiguration();
        AppUtils::logger_ptr logger();
        std::string resultName() const;
        std::string trajectoriesResultName() const;
        std::string confFileName() const;
        int numberOfThreads() const;

    private:
        simConf_ptr simulationConfiguration_;
        AppUtils::logger_ptr logger_;
        std::string simResultName_;
        std::string confFileName_;
        int numberOfThreads_= 1;
    };
}




#endif //IDSIMF_APPUTILS_COMMANDLINEPARSER_HPP
