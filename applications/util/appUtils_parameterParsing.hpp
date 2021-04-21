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
        bool isParameter(const std::string& keyName);
        int intConfParameter(const std::string& jsonName);
        std::vector<int> intVectorParameter(const std::string& jsonName);
        double doubleParameter(const std::string& jsonName);
        std::vector<double> doubleVectorParameter(const std::string& jsonName, double multiplicator = 1.0);
        Core::Vector vector3dConfParameter(const std::string& jsonName);
        std::array<std::array<double,2>,3> double3dBox(const std::string& jsonName);
        std::string stringConfParameter(const std::string& jsonName);
        std::vector<std::string> stringVectorConfParameter(const std::string& jsonName);
        bool boolConfParameter(const std::string& jsonName);

        std::unique_ptr<ParticleSimulation::InterpolatedField> readInterpolatedField(
                const std::string& jsonName);

        std::string pathRelativeToConfFile(const std::string& pathStr);
        std::string pathRelativeToConfBasePath(const std::string& pathStr);

    private:
        Json::Value readConfigurationJson_(const std::string& confFileName);

        Json::Value confRoot_;
        std::filesystem::path confFilePath_;
        std::filesystem::path confFileBasePath_;
    };
}
#endif //IONSIMULATION_CPP_PARAMETERPARSING_HPP
