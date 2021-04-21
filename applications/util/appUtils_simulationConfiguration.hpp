#ifndef IONSIMULATION_CPP_PARAMETERPARSING_HPP
#define IONSIMULATION_CPP_PARAMETERPARSING_HPP

#include "json.h"
#include "PSim_interpolatedField.hpp"
#include <filesystem>
#include <vector>
#include <fstream>
#include <iostream>

namespace AppUtils{
    class SimulationConfiguration{
    public:
        SimulationConfiguration(const std::string& confFileName);
        bool isParameter(const std::string& keyName) const;
        int intParameter(const std::string& jsonName) const;
        std::vector<int> intVectorParameter(const std::string& jsonName) const;
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
    };
}
#endif //IONSIMULATION_CPP_PARAMETERPARSING_HPP
