#ifndef IONSIMULATION_CPP_PARAMETERPARSING_HPP
#define IONSIMULATION_CPP_PARAMETERPARSING_HPP

#include "json.h"
#include "PSim_interpolatedField.hpp"
#include <filesystem>
#include <vector>
#include <fstream>
#include <iostream>

Json::Value readConfigurationJson(std::string confFileName);
bool isConfFileKey(std::string keyName, Json::Value& confRoot);
int intConfParameter(std::string jsonName, Json::Value& confRoot);
std::vector<int> intVectorConfParameter(std::string jsonName, Json::Value& confRoot);
double doubleConfParameter(std::string jsonName, Json::Value& confRoot);
std::vector<double> doubleVectorConfParameter(std::string jsonName, Json::Value& confRoot, double multiplicator = 1.0);
std::array<std::array<double,2>,3> double3dBox(std::string jsonName, Json::Value& confRoot);
std::string stringConfParameter(std::string jsonName, Json::Value& confRoot);
std::vector<std::string> stringVectorConfParameter(std::string jsonName, Json::Value& confRoot);

bool boolConfParameter(std::string jsonName, Json::Value& confRoot);

std::unique_ptr<ParticleSimulation::InterpolatedField> readInterpolatedField(
        std::string confBasePathStr,
        std::string jsonName,
        Json::Value& confRoot);

std::string pathRelativeToConfFile(std::string confFilePathStr, std::string pathStr);
std::string pathRelativeToConfBasePath(std::string confBasePath, std::string pathStr);

std::string confFileBasePath(std::string confFilePathStr);


#endif //IONSIMULATION_CPP_PARAMETERPARSING_HPP
