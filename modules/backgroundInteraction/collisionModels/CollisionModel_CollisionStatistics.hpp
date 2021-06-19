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

 ------------
 CollisionModel_CollisionStatistics.hpp

 A collision statistics: Describes the normalized statistics of
 diffusive particle migration due to collisions in terms of
 inverse cumulative density functions (ICDFs).

 ****************************/

#ifndef Collision_CollisionStatistics_hpp
#define Collision_CollisionStatistics_hpp

#include <stdexcept>
#include <vector>
#include <string>

namespace CollisionModel{
    /**
     * Individual exception class for problems in collision statistics files
    */
    class CollisionStatisticsFileException : public std::runtime_error {
    public:
        explicit CollisionStatisticsFileException (const std::string msg): std::runtime_error(msg) {}
    };

    class CollisionStatistics {
        public:

            CollisionStatistics();
            explicit CollisionStatistics(std::string statisticsFilename);
            int getNDist() const;
            int getNDistPoints() const;
            int getNCollisions() const;
            std::vector<double> getMassRatios() const;
            double getLogMassRatio(std::size_t icdfIndex) const;
            double getLogMassRatioDistance(std::size_t index) const;
            std::vector<std::vector<double>> getICDFs() const;
            std::size_t findUpperDistIndex(const double logMassRatio) const;


        private:
            int nDist_ = 0; ///< Number of ICDFs in represented in distribution data
            int nDistPoints_ = 0; ///< Number of data points per ICDF.
            int nDistCollisions_ = 0; ///< Number of collisions represented in distribution data
            std::vector<double> massRatios_; ///< Mass ratios in the individual ICDFs
            std::vector<double> logMassRatios_; ///< Logarithm of the mass ratios
            std::vector<double> logMassRatioDists_; ///< Logarithmic distances between the mass ratios
            std::vector<std::vector<double>> icdfs_; ///< The ICDF data in the collision statistic

            void initCollisionModelFromStatisticsFile_(std::string statisticsFilename);
    };
}

#endif //Collision_CollisionStatistic_hpp
