#ifndef IONSIMULATION_CPP_PARAMETERPARSING_HPP
#define IONSIMULATION_CPP_PARAMETERPARSING_HPP

#include "json.h"
#include "PSim_interpolatedField.hpp"
#include <filesystem>
#include <vector>
#include <fstream>
#include <iostream>

Json::Value readConfigurationJson(const std::string& confFileName);
bool isConfFileKey(const std::string& keyName, const Json::Value& confRoot);
int intConfParameter(const std::string& jsonName, const Json::Value& confRoot);
std::vector<int> intVectorConfParameter(const std::string& jsonName, const Json::Value& confRoot);
double doubleConfParameter(const std::string& jsonName, const Json::Value& confRoot);
std::vector<double> doubleVectorConfParameter(const std::string& jsonName, const Json::Value& confRoot, double multiplicator = 1.0);
Core::Vector vector3dConfParameter(const std::string& jsonName, const Json::Value& confRoot);
std::array<std::array<double,2>,3> double3dBox(const std::string& jsonName, const Json::Value& confRoot);
std::string stringConfParameter(const std::string& jsonName, const Json::Value& confRoot);
std::vector<std::string> stringVectorConfParameter(const std::string& jsonName, const Json::Value& confRoot);

bool boolConfParameter(const std::string& jsonName, const Json::Value& confRoot);

std::unique_ptr<ParticleSimulation::InterpolatedField> readInterpolatedField(
        const std::string& confBasePathStr,
        const std::string& jsonName,
        const Json::Value& confRoot);

std::string pathRelativeToConfFile(const std::string& confFilePathStr, const std::string& pathStr);
std::string pathRelativeToConfBasePath(const std::string& confBasePath, const std::string& pathStr);

std::string confFileBasePath(const std::string& confFilePathStr);


#endif //IONSIMULATION_CPP_PARAMETERPARSING_HPP
