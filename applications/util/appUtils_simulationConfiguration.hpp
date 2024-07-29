/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2024 - Physical and Theoretical Chemistry /
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
 appUtils_simulationConfiguration.hpp

 Simulation configuration class:
 A simulation configuration defines the configuration of a simulation run from the simulation input / configuration file
 ****************************/
#ifndef IONSIMULATION_CPP_PARAMETERPARSING_HPP
#define IONSIMULATION_CPP_PARAMETERPARSING_HPP

#include "json.h"
#include "spdlog/spdlog.h"
#include "PSim_interpolatedField.hpp"
#include "appUtils_logging.hpp"
#include <filesystem>
#include <vector>
#include <fstream>
#include <iostream>
#include <memory>

namespace AppUtils{

    enum IntegratorMode {VERLET, PARALLEL_VERLET, FMM3D_VERLET, EXAFMM_VERLET, FULL_SUM_VERLET, PARALLEL_RUNGE_KUTTA4, FULL_SUM_RUNGE_KUTTA4};

    class SimulationConfiguration{
    public:
        SimulationConfiguration(const std::string& confFileName);
        SimulationConfiguration(const std::string& confFileName, std::shared_ptr<spdlog::logger> logger);

        bool isParameter(const std::string& keyName) const;
        bool isVectorParameter(const std::string& keyName) const;

        int intParameter(const std::string& jsonName) const;
        unsigned int unsignedIntParameter(const std::string& jsonName) const;
        std::vector<int> intVectorParameter(const std::string& jsonName) const;
        std::vector<unsigned int> unsignedIntVectorParameter(const std::string& jsonName) const;

        double doubleParameter(const std::string& jsonName) const;
        std::vector<double> doubleVectorParameter(const std::string& jsonName, double multiplicator = 1.0) const;
        Core::Vector vector3dParameter(const std::string& jsonName) const;
        std::array<std::array<double,2>,3> double3dBox(const std::string& jsonName) const;

        std::string stringParameter(const std::string& jsonName) const;
        std::vector<std::string> stringVectorParameter(const std::string& jsonName) const;

        bool boolParameter(const std::string& jsonName) const;

        AppUtils::IntegratorMode integratorMode() const;
        std::unique_ptr<ParticleSimulation::InterpolatedField> readInterpolatedField(
                const std::string& jsonName) const;

        std::string pathRelativeToConfFile(const std::string& pathStr) const;
        std::string pathRelativeToConfBasePath(const std::string& pathStr) const;
        std::string confBasePath() const;

    private:
        Json::Value readConfigurationJson_(const std::string& confFileName);

        Json::Value confRoot_;
        std::filesystem::path confFilePath_;
        std::filesystem::path confFileBasePath_;
        AppUtils::logger_ptr logger_ = nullptr;
    };

    using simConf_ptr = std::shared_ptr<AppUtils::SimulationConfiguration>;
}
#endif //IONSIMULATION_CPP_PARAMETERPARSING_HPP
