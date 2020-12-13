#include "appUtils_parameterParsing.hpp"
#include "Core_vector.hpp"

Json::Value readConfigurationJson(const std::string& confFileName){
    std::cout << confFileName<<std::endl;
    std::filesystem::path confFilePath(confFileName);
    if (!std::filesystem::exists(confFilePath)){
        throw std::invalid_argument("Configuration file not existing: " + confFilePath.string());
    }

    std::ifstream confFile;
    confFile.open(confFilePath);
    Json::Value confRoot;
    confFile >> confRoot;
    confFile.close();
    std::cout<<confRoot<<std::endl;

    return confRoot;
}

bool isConfFileKey(const std::string& keyName, const Json::Value& confRoot) {
    return confRoot.isMember(keyName);
}

int intConfParameter(const std::string& jsonName, const Json::Value& confRoot){
    if (confRoot.isMember(jsonName)) {
        int result = confRoot.get(jsonName, 0).asInt();
        std::cout << jsonName << ":" << result << std::endl;
        return(result);
    } else {
        throw std::invalid_argument("missing configuration value: " + jsonName);
    }
}

std::vector<int> intVectorConfParameter(const std::string& jsonName, const Json::Value& confRoot){
    std::vector<int> result;
    if (confRoot.isMember(jsonName)) {
        Json::Value jsonNode = confRoot.get(jsonName,0);
        for (int i=0; i<jsonNode.size(); i++){
            result.push_back(jsonNode.get(i,0).asInt());
        }
    } else {
        throw std::invalid_argument("missing configuration value: " + jsonName);
    }
    return (result);
}

double doubleConfParameter(const std::string& jsonName, const Json::Value& confRoot){
    if (confRoot.isMember(jsonName) == true) {
        double result = confRoot.get(jsonName, 0).asDouble();
        std::cout << jsonName << ":" << result << std::endl;
        return(result);
    } else {
        throw std::invalid_argument("missing configuration value: " + jsonName);
    }
}

std::vector<double> doubleVectorConfParameter(const std::string& jsonName, const Json::Value& confRoot, double multiplicator){
    std::vector<double> result;
    if (confRoot.isMember(jsonName) == true) {
        Json::Value jsonNode = confRoot.get(jsonName,0);
        for (int i=0; i<jsonNode.size(); i++){
            result.push_back(jsonNode.get(i,0).asDouble()*multiplicator);
        }
    } else {
        throw std::invalid_argument("missing configuration value: " + jsonName);
    }
    return (result);
}

Core::Vector vector3dConfParameter(const std::string& jsonName, const Json::Value& confRoot){
    std::vector<double> vectorRaw = doubleVectorConfParameter(jsonName, confRoot);
    return {vectorRaw[0], vectorRaw[1], vectorRaw[2]};
}

std::array<std::array<double,2>,3> double3dBox(const std::string& jsonName, const Json::Value& confRoot){
    std::array<std::array<double, 2>, 3> result{{{0, 0}, {0, 0}, {0, 0}}};
    if (confRoot.isMember(jsonName) == true) {
        Json::Value jsonNode = confRoot.get(jsonName,0);

        if (jsonNode.size() != 3){
            throw std::invalid_argument("3d box definition has not 3 dimensions");
        }

        for (int i=0; i<3; ++i){
            Json::Value dimNode = jsonNode.get(i,0);
            if (dimNode.size() != 2){
                throw std::invalid_argument("Number of elements is not 2 in a dimension of a 3d box definition");
            }
            result[i][0] = dimNode.get(int(0) ,0).asDouble();
            result[i][1] = dimNode.get(int(1) ,0).asDouble();
        }
    } else {
        throw std::invalid_argument("missing configuration value: " + jsonName);
    }
    return (result);
}

std::string stringConfParameter(const std::string& jsonName, const Json::Value& confRoot){
    if (confRoot.isMember(jsonName) == true) {
        std::string result = confRoot.get(jsonName, 0).asString();
        std::cout << jsonName << ":" << result << std::endl;
        return(result);
    } else {
        throw std::invalid_argument("missing configuration value: " + jsonName);
    }
}

std::vector<std::string> stringVectorConfParameter(const std::string& jsonName, const Json::Value& confRoot){
    std::vector<std::string> result;
    if (confRoot.isMember(jsonName) == true) {
        Json::Value jsonNode = confRoot.get(jsonName,0);
        for (int i=0; i<jsonNode.size(); i++){
            result.push_back(jsonNode.get(i,0).asString());
        }
    } else {
        throw std::invalid_argument("missing configuration value: " + jsonName);
    }
    return (result);
}

bool boolConfParameter(const std::string& jsonName, const Json::Value& confRoot){
    if (confRoot.isMember(jsonName) == true) {
        std::string confString = confRoot.get(jsonName, 0).asString();
        std::cout << jsonName << ":" << confString << std::endl;
        if (confString == "true") {
            return true;
        }
        else if (confString == "false") {
            return false;
        }
        else {
            throw std::invalid_argument("invalid boolean option: " + jsonName);
        }
    } else {
        throw std::invalid_argument("missing configuration value: " + jsonName);
    }
}

std::unique_ptr<ParticleSimulation::InterpolatedField> readInterpolatedField(
        const std::string& confBasePathStr,
        const std::string& jsonName,
        const Json::Value& confRoot){
    if (confRoot.isMember(jsonName)==true){
        std::string fieldFileName = confRoot.get(jsonName,0).asString();
        std::filesystem::path fieldPath(std::filesystem::path(confBasePathStr) / std::filesystem::path(fieldFileName));
        std::cout << "Reading field "<<fieldPath <<" \n";
        auto fieldPtr = std::make_unique<ParticleSimulation::InterpolatedField>(fieldPath);
        return move(fieldPtr);
    }else{
        throw std::invalid_argument("missing configuration value: " + jsonName);
    }
}

std::string pathRelativeToConfFile(const std::string& confFilePathStr, const std::string& pathStr){
    return std::filesystem::path(confFilePathStr).parent_path() / std::filesystem::path(pathStr);
}

std::string pathRelativeToConfBasePath(const std::string& confBasePath, const std::string& pathStr){
    return std::filesystem::path(confBasePath) / std::filesystem::path(pathStr);
}

std::string confFileBasePath(const std::string& confFilePathStr){
    return std::filesystem::path(confFilePathStr).parent_path();
}