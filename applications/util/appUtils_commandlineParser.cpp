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

AppUtils::CommandlineParser::CommandlineParser(
        int argc, const char** argv, std::string appName, std::string appDescription,  bool multithreaded) {
    CLI::App app{appDescription, appName};

    app.add_option("--conf,-c,run_config",confFileName_, "Configuration file")->required();
    app.add_option("--result,-r,result",projectName_, "Result name")->required();

    if (multithreaded){
        app.add_option("--n_threads,-n", numberOfThreads_, "number of parallel threads")->required();
    }

    int retCode = parse_(app, argc, argv);
    if (retCode != 0){
        throw std::runtime_error("Wrong commandline");
    }

    if (multithreaded){
        omp_set_num_threads(numberOfThreads_);
    }
    logger_ = AppUtils::createLogger(projectName_ + ".log");
    simulationConfiguration_ = std::make_shared<SimulationConfiguration>(confFileName_, logger_);
}

int AppUtils::CommandlineParser::parse_(CLI::App& app, int argc, const char * argv[]) {
    CLI11_PARSE(app, argc, argv);
    return 0;
}

std::shared_ptr<AppUtils::SimulationConfiguration> AppUtils::CommandlineParser::simulationConfiguration() {
    return simulationConfiguration_;
}

std::shared_ptr<spdlog::logger> AppUtils::CommandlineParser::logger() {
    return logger_;
}

std::string AppUtils::CommandlineParser::confFileName() {
    return confFileName_;
}

std::string AppUtils::CommandlineParser::projectName() {
    return projectName_;
}

int AppUtils::CommandlineParser::numberOfThreads() {
    return  numberOfThreads_;
}