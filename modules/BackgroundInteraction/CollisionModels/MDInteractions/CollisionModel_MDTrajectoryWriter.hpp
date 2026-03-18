/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2025 - Physical and Theoretical Chemistry /
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
 CollisionModel_MDTrajectoryWriter.hpp

 A writer for MD collision trajectories

 ****************************/
#ifndef COLLISIONMODEL_MDTRAJECTORYWRITER_HPP
#define COLLISIONMODEL_MDTRAJECTORYWRITER_HPP

#include "Core_vector.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>

namespace CollisionModel{
    class MDTrajectoryWriter {
        public:
            explicit MDTrajectoryWriter(const std::string& trajectoryFileName, double minimalSampleInterval = 0.0);

            void writeTrajectorySample(double time, double dt,
                                       Core::Vector positionBgMolecule,
                                       Core::Vector velocityBgMolecule, Core::Vector positionMolecule, std::vector<Core::Vector> forceMolecules,
                                       double distance);

            void writeTrajectoryDelimiter();


        private:
            std::unique_ptr<std::ofstream> trajectoryOutputStream_;
            double minimalSampleInterval_;
            double nextTrajectorySampleTime_ = 0.0;
    };
}


#endif //COLLISIONMODEL_MDTRAJECTORYWRITER_HPP
