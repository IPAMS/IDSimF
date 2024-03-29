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

#include "CollisionModel_CollisionStatistics.hpp"
#include "CollisionStatistic_default.hpp"
#include "CollisionModel_MathFunctions.hpp"
#include <algorithm>
#include <regex>

/**
 * Default constructor: Inits collision statistics with a default dataset
 */
CollisionModel::CollisionStatistics::CollisionStatistics():
        nDist_(CollisionModel::N_DIST_DEFAULT),
        nDistPoints_(CollisionModel::N_DIST_DEFAULT_POINTS),
        nDistCollisions_(CollisionModel::N_DIST_DEFAULT_COLLISIONS),
        massRatios_(CollisionModel::MASSRATIO_DEFAULT),
        logMassRatios_(CollisionModel::LOG_MASSRATIO_DEFAULT),
        logMassRatioDists_(CollisionModel::LOG_MASSRATIO_DISTANCES_DEFAULT),
        icdfs_(CollisionModel::ICDF_DEFAULT)
{}

/**
 * Constructor: Creates collision statistics from a collision statistics file
 * (e.g. generated by the CollisionStatistics.jl)
 *
 * @param statisticsFilename the name / path to a collision statistics parameter file
 */
CollisionModel::CollisionStatistics::CollisionStatistics(const std::string &statisticsFilename) {
    initCollisionModelFromStatisticsFile_(statisticsFilename);
    std::transform(massRatios_.begin(), massRatios_.end(), std::back_inserter(logMassRatios_),
            [](double mr) -> double { return log10(mr); });
    for (size_t i=1; i< logMassRatios_.size(); ++i){
        logMassRatioDists_.emplace_back(logMassRatios_[i]-logMassRatios_[i-1]);
    }
}

/**
 * Inits the collision statistics from a collision statistics file
 * as exported from CollisionStatistics.jl
 *
 * @param statisticsFilename the name / path to a collision statistics parameter file
 */
void CollisionModel::CollisionStatistics::initCollisionModelFromStatisticsFile_(const std::string &statisticsFilename) {
    std::ifstream instream;
    instream.open(statisticsFilename);

    if (instream.good()){

        std::regex patternNumericParameter(R"(;\s*([\w_ ]+)\s*=\s*([\d]+))");
        std::regex patternNumber(R"(^\d+.?\d*\r?$)");

        std::smatch matches;

        for (std::string line; std::getline(instream, line); ) {
            if (std::regex_match(line, patternNumber)){
                icdfs_.back().push_back(std::strtof(line.c_str(), nullptr));
            }
            else if (std::regex_search(line, matches, patternNumericParameter)) {
                std::string parameterName = matches[1];
                std::string parameterValue = matches[2];

                if (parameterName == "ICDF_massratio"){
                    massRatios_.push_back(std::strtof(parameterValue.c_str(), nullptr));
                    icdfs_.emplace_back(std::vector<double>(0));
                }
                else if (parameterName == "n_collisions"){
                    nDistCollisions_ = static_cast<int>(std::stol(parameterValue));
                }
                else if (parameterName == "n_dist_points"){
                    nDistPoints_ = static_cast<int>(std::stol(parameterValue));
                }
                else if (parameterName == "n_statistics"){
                    nDist_ = static_cast<int>(std::stol(parameterValue));
                }
            }
        }
    }
    else {
        throw(CollisionModel::CollisionStatisticsFileException("file not found"));
    }
}

/**
 * Gets the index of the collision distribution with a logMassRatio above the given logMassRatio
 * @param logMassRatio the logarithmic mass ratio to find the distribution index for
 * @return index of the collision distribution with a logarithmic mass ratio above logMassRatio
 */
std::size_t CollisionModel::CollisionStatistics::findUpperDistIndex(const double logMassRatio) const {
    auto upper = std::upper_bound(logMassRatios_.begin(), logMassRatios_.end(), logMassRatio);
    if (upper == logMassRatios_.begin()){
        return 0;
    }
    else {
        return ( static_cast<std::size_t>((--upper) - logMassRatios_.begin()));
    }
}


/**
 * Gets the number of collision statistics (individual collision ICDFs)
 *
 * @return the number of individual collision ICDFs in the collision statistics
 */
int CollisionModel::CollisionStatistics::getNDist() const {
    return nDist_;
}

/**
 * Gets the number of points / steps in the probabilty dimension of the
 * ICDFs. The probability interval [0..1] is splitted in nDistPoints equidistant
 * steps.
 *
 * @return Number of points / steps in the ICDFs.
 */
int CollisionModel::CollisionStatistics::getNDistPoints() const {
    return nDistPoints_;
}

/**
 * Gets the number of collisions which was performed to simulate / calculate the
 * ICDFs.
 *
 * @return Number of collisions in the Monte-Carlo simulations leading to the ICDFs.
 */
int CollisionModel::CollisionStatistics::getNCollisions() const {
    return nDistCollisions_;
}

/**
 * Returns the mass ratios of the ICDFs.
 *
 * @return A vector with the mass ratios used in the Monte-Carlo simulations leading to
 * the ICDFs.
 */
std::vector<double> CollisionModel::CollisionStatistics::getMassRatios() const {
    return massRatios_;
}

/**
 * Returns the logarithm of the mass ratio of an ICDF
 *
 * @param icdfIndex the index of the ICDF
 * @return log(mass ratio) of the ICDF with index icdfIndex
 */
double CollisionModel::CollisionStatistics::getLogMassRatio(size_t icdfIndex) const {
    return logMassRatios_[icdfIndex];
}

/**
 * Returns the difference of the logarithmic mass ratios between
 * ICDF with indices index and index+1
 * @param index the index of the mass ratio difference
 * @return difference of logarithmic mass ratios of ICDF with indices index and index+1
 */
double CollisionModel::CollisionStatistics::getLogMassRatioDistance(size_t index) const {
    return logMassRatioDists_[index];
}

/**
 * Gets the ICDFs from the collision statistics.
 * The ICDFs are a vector of double vectors, one per mass ratio, with nDistPoints steps each.
 *
 * @return The ICDFs from the collision statistics
 */
std::vector<std::vector<double>> CollisionModel::CollisionStatistics::getICDFs() const {
    return icdfs_;
}