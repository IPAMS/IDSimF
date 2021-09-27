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
 ****************************/

#include "appUtils_commandlineParser.hpp"
#include "appUtils_logging.hpp"
#include <omp.h>

/**
 * Parses the commandline options / switches and constructs an object with the simulation configuration,
 * simulation result name and additional options specified by the commandline.
 *
 * @param argc Number of commandline arguments
 * @param argv Commandline argument list
 * @param appName Name of the simulation application
 * @param appDescription A brief description of the simualtion application
 * @param multithreaded If true the commandline arguments are parsed for a multithreaded application.
 */
AppUtils::CommandlineParser::CommandlineParser(
        int argc, const char** argv, std::string appName, std::string appDescription,  bool multithreaded) {
    CLI::App app{appDescription, appName};

    app.add_option("--conf,-c,run_config",confFileName_, "Configuration file")->required();
    app.add_option("--result,-r,result",simResultName_, "Result name")->required();

    if (multithreaded){
        app.add_option("--n_threads,-n", numberOfThreads_, "number of parallel threads")->required();
    }

    // Try to parse and raise a message to the main app if the parsing fails for some reason:
    try {
        (app).parse((argc), (argv));
    } catch(const CLI::ParseError &e) {
        int returnCode = (app).exit(e);
        throw AppUtils::TerminatedWhileCommandlineParsing(returnCode, "app terminated while commandline parsing");
    }

    if (multithreaded){
        omp_set_num_threads(numberOfThreads_);
    }
    logger_ = AppUtils::createLogger(simResultName_ + ".log");
    simulationConfiguration_ = std::make_shared<SimulationConfiguration>(confFileName_, logger_);
}

/**
 * Returns the (parsed) simulation configuration specified by the simulation configuration file name in the commandline
 */
std::shared_ptr<AppUtils::SimulationConfiguration> AppUtils::CommandlineParser::simulationConfiguration() {
    return simulationConfiguration_;
}

/**
 * Returns the logger object for the app, the log file name is constructed from the result base name in the commandline
 */
std::shared_ptr<spdlog::logger> AppUtils::CommandlineParser::logger() {
    return logger_;
}

/**
 * Returns the name of the configuration file name specified in the commandline
 */
std::string AppUtils::CommandlineParser::confFileName() {
    return confFileName_;
}

/**
 * Returns the result base name specified in the commandline
 */
std::string AppUtils::CommandlineParser::resultName() {
    return simResultName_;
}

/**
 * Returns the number of parallel threads specified by the commandline. For a single threaded application, the number
 * of threads is always 1
 */
int AppUtils::CommandlineParser::numberOfThreads() {
    return  numberOfThreads_;
}