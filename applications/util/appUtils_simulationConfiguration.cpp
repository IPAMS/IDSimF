#include "appUtils_simulationConfiguration.hpp"
#include "appUtils_inputFileUtilities.hpp"
#include "Core_vector.hpp"

AppUtils::SimulationConfiguration::SimulationConfiguration(const std::string& confFileName) {

    confRoot_ = readConfigurationJson_(confFileName);
}

AppUtils::SimulationConfiguration::SimulationConfiguration(const std::string& confFileName,
                                                           std::shared_ptr<spdlog::logger> logger):
    logger_(logger)
{
    confRoot_ = readConfigurationJson_(confFileName);
}

bool AppUtils::SimulationConfiguration::isParameter(const std::string& keyName) const {
    return confRoot_.isMember(keyName);
}

bool AppUtils::SimulationConfiguration::isVectorParameter(const std::string& keyName) const {
    Json::Value value = confRoot_.get(keyName, 0);
    return value.isArray();
}

int AppUtils::SimulationConfiguration::intParameter(const std::string& jsonName) const {
    if (isParameter(jsonName)) {
        int result = confRoot_.get(jsonName, 0).asInt();
        if (logger_){
            logger_->info("{}:{}",jsonName, result);
        }
        return(result);
    } else {
        throw std::invalid_argument("missing configuration value: " + jsonName);
    }
}

unsigned int AppUtils::SimulationConfiguration::unsignedIntParameter(const std::string& jsonName) const {
    if (isParameter(jsonName)) {
        unsigned int result = confRoot_.get(jsonName, 0).asUInt();
        if (logger_){
            logger_->info("{}:{}",jsonName, result);
        }
        return(result);
    } else {
        throw std::invalid_argument("missing configuration value: " + jsonName);
    }
}

std::vector<int> AppUtils::SimulationConfiguration::intVectorParameter(const std::string& jsonName) const {
    std::vector<int> result;
    if (isParameter(jsonName)) {
        Json::Value jsonNode = confRoot_.get(jsonName,0);
        auto nJsonNodes = static_cast<unsigned int>(jsonNode.size());
        for (unsigned int i=0; i<nJsonNodes; i++){
            result.push_back(jsonNode.get(i,0).asInt());
        }
    } else {
        throw std::invalid_argument("missing configuration value: " + jsonName);
    }
    return (result);
}

std::vector<unsigned int> AppUtils::SimulationConfiguration::unsignedIntVectorParameter(const std::string& jsonName) const {
    std::vector<unsigned int> result;
    if (isParameter(jsonName)) {
        Json::Value jsonNode = confRoot_.get(jsonName,0);
        auto nJsonNodes = static_cast<unsigned int>(jsonNode.size());

        for (unsigned int i=0; i<nJsonNodes; ++i){
            result.push_back(jsonNode.get(i,0).asUInt());
        }
    } else {
        throw std::invalid_argument("missing configuration value: " + jsonName);
    }
    return (result);
}


double AppUtils::SimulationConfiguration::doubleParameter(const std::string& jsonName) const {
    if (isParameter(jsonName)) {
        double result = confRoot_.get(jsonName, 0).asDouble();
        if (logger_){
            logger_->info("{}:{}",jsonName, result);
        }
        return(result);
    } else {
        throw std::invalid_argument("missing configuration value: " + jsonName);
    }
}

std::vector<double> AppUtils::SimulationConfiguration::doubleVectorParameter(const std::string& jsonName,
                                                                             double multiplicator) const {
    std::vector<double> result;
    if (isParameter(jsonName)) {
        Json::Value jsonNode = confRoot_.get(jsonName,0);
        auto nJsonNodes = static_cast<unsigned int>(jsonNode.size());
        for (unsigned int i=0; i<nJsonNodes; ++i){
            result.push_back(jsonNode.get(i,0).asDouble()*multiplicator);
        }
    } else {
        throw std::invalid_argument("missing configuration value: " + jsonName);
    }
    return (result);
}

Core::Vector AppUtils::SimulationConfiguration::vector3dParameter(const std::string& jsonName) const {
    std::vector<double> vectorRaw = doubleVectorParameter(jsonName);
    Core::Vector result{vectorRaw[0], vectorRaw[1], vectorRaw[2]};
    if (logger_){
        logger_->info("{}:{} {} {}",jsonName, result.x(), result.y(), result.z());
    }
    return result;
}


std::array<std::array<double, 2>, 3> AppUtils::SimulationConfiguration::double3dBox(const std::string& jsonName) const {
    std::array<std::array<double, 2>, 3> result{{{0, 0}, {0, 0}, {0, 0}}};
    if (isParameter(jsonName)) {
        Json::Value jsonNode = confRoot_.get(jsonName,0);

        if (jsonNode.size() != 3){
            throw std::invalid_argument("3d box definition has not 3 dimensions");
        }

        for (unsigned int i=0; i<3; ++i){
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

std::string AppUtils::SimulationConfiguration::stringParameter(const std::string& jsonName) const {
    if (isParameter(jsonName)) {
        std::string result = confRoot_.get(jsonName, 0).asString();
        if (logger_){
            logger_->info("{}:{}",jsonName, result);
        }
        return(result);
    } else {
        throw std::invalid_argument("missing configuration value: " + jsonName);
    }
}

std::vector<std::string> AppUtils::SimulationConfiguration::stringVectorParameter(const std::string& jsonName) const {
    std::vector<std::string> result;
    if (isParameter(jsonName)) {
        Json::Value jsonNode = confRoot_.get(jsonName,0);
        auto nJsonNodes = static_cast<unsigned int>(jsonNode.size());
        for (unsigned int i=0; i<nJsonNodes; ++i){
            result.push_back(jsonNode.get(i,0).asString());
        }
    } else {
        throw std::invalid_argument("missing configuration value: " + jsonName);
    }
    return (result);
}

bool AppUtils::SimulationConfiguration::boolParameter(const std::string& jsonName) const {
    if (isParameter(jsonName)) {
        std::string confString = confRoot_.get(jsonName, 0).asString();
        if (logger_){
            logger_->info("{}:{}",jsonName, confString);
        }
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

AppUtils::IntegratorMode AppUtils::SimulationConfiguration::integratorMode() const {
    std::string integratorMode_str = this->stringParameter("integrator_mode");
    IntegratorMode integratorMode;
    if (integratorMode_str=="verlet") {
        integratorMode = VERLET;
    }
    else if (integratorMode_str=="parallel_verlet") {
        integratorMode = PARALLEL_VERLET;
    }
    else if (integratorMode_str=="RK4") {
        integratorMode = PARALLEL_RUNGE_KUTTA4;
    }
    else if (integratorMode_str=="full_sum_RK4") {
        integratorMode = FULL_SUM_RUNGE_KUTTA4;
    }
    else if (integratorMode_str=="full_sum_verlet") {
        integratorMode = FULL_SUM_VERLET;
    }
#ifdef WITH_FMM_3d
    else if (integratorMode_str=="FMM3D_verlet") {
        integratorMode = FMM3D_VERLET;
    }
#endif
#ifdef WITH_EXAFMMT
    else if (integratorMode_str=="ExaFMM_verlet") {
        integratorMode = EXAFMM_VERLET;
    }
#endif
    else {
        throw std::invalid_argument("wrong configuration value: integrator mode");
    }
    return integratorMode;
}

std::unique_ptr<ParticleSimulation::InterpolatedField> AppUtils::SimulationConfiguration::readInterpolatedField(
        const std::string& jsonName) const {
        if (isParameter(jsonName)){
            std::string fieldFileName = confRoot_.get(jsonName,0).asString();
            std::filesystem::path fieldPath(confFileBasePath_ / std::filesystem::path(fieldFileName));
            if (logger_){
                logger_->info("Reading field {}", fieldPath.string());
            }

            return std::make_unique<ParticleSimulation::InterpolatedField>(fieldPath);
        }else{
            throw std::invalid_argument("missing configuration value: " + jsonName);
        }
}

std::vector<std::unique_ptr<ParticleSimulation::SimionPotentialArray> > AppUtils::SimulationConfiguration::readPotentialArrays(
    const std::string& jsonName, double paSpatialScale, bool fastAdjustPA) const {
    if (isParameter(jsonName)){
        std::vector<std::string> potentialArraysNames = stringVectorParameter(jsonName);
        std::vector<std::unique_ptr<ParticleSimulation::SimionPotentialArray>> potentialArrays =
            readPotentialArrayFiles(potentialArraysNames, confBasePath(), paSpatialScale, fastAdjustPA);

        return potentialArrays;
    }else{
        throw std::invalid_argument("missing configuration value: " + jsonName);
    }
}

std::string AppUtils::SimulationConfiguration::pathRelativeToConfFile(const std::string& pathStr) const {
    return confFilePath_.parent_path() / std::filesystem::path(pathStr);
}

std::string AppUtils::SimulationConfiguration::pathRelativeToConfBasePath(const std::string& pathStr) const {
    return confFileBasePath_ / std::filesystem::path(pathStr);
}

std::string AppUtils::SimulationConfiguration::confBasePath() const {
    return confFileBasePath_;
}

Json::Value AppUtils::SimulationConfiguration::readConfigurationJson_(const std::string& confFileName) {
    if (logger_){
        logger_->info("Reading configuration file:{}", confFileName);
    }
    confFilePath_ = std::filesystem::path(confFileName);
    if (!std::filesystem::exists(confFilePath_)){
        throw std::invalid_argument("Configuration file not existing: " + confFilePath_.string());
    }
    confFileBasePath_ = confFilePath_.parent_path();

    std::ifstream confFile;
    confFile.open(confFilePath_);
    Json::Value confRoot;
    confFile >> confRoot;
    confFile.close();
    return confRoot;
}