#ifndef IONSIMULATION_CPP_PARAMETERPARSING_HPP
#define IONSIMULATION_CPP_PARAMETERPARSING_HPP

#include "json.h"
#include "spdlog/spdlog.h"
#include "PSim_interpolatedField.hpp"
#include <filesystem>
#include <vector>
#include <fstream>
#include <iostream>
#include <memory>

namespace AppUtils{
    class SimulationConfiguration{
    public:
        SimulationConfiguration(const std::string& confFileName);
        SimulationConfiguration(const std::string& confFileName, std::shared_ptr<spdlog::logger> logger);
        bool isParameter(const std::string& keyName) const;
        int intParameter(const std::string& jsonName) const;
        std::vector<int> intVectorParameter(const std::string& jsonName) const;
        std::vector<unsigned int> unsignedIntVectorParameter(const std::string& jsonName) const;
        double doubleParameter(const std::string& jsonName) const;
        std::vector<double> doubleVectorParameter(const std::string& jsonName, double multiplicator = 1.0) const;
        Core::Vector vector3dParameter(const std::string& jsonName) const;
        std::array<std::array<double,2>,3> double3dBox(const std::string& jsonName) const;
        std::string stringParameter(const std::string& jsonName) const;
        std::vector<std::string> stringVectorParameter(const std::string& jsonName) const;
        bool boolParameter(const std::string& jsonName) const;

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
        std::shared_ptr<spdlog::logger> logger_ = nullptr;
    };
}
#endif //IONSIMULATION_CPP_PARAMETERPARSING_HPP
