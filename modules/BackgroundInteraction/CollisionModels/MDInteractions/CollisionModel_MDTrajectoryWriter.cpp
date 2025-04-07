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
 ****************************/
#include "CollisionModel_MDTrajectoryWriter.hpp"


CollisionModel::MDTrajectoryWriter::MDTrajectoryWriter(const std::string& trajectoryFileName, double minimalSampleInterval) {
    trajectoryOutputStream_ = std::make_unique<std::ofstream>();
    trajectoryOutputStream_->open(trajectoryFileName, std::ofstream::app);

    if (trajectoryOutputStream_->good()){
        minimalSampleInterval_ = minimalSampleInterval;
    }
    else{
        throw (std::runtime_error("Trajectory Output Stream failed to open"));
    }
}

void CollisionModel::MDTrajectoryWriter::writeTrajectorySample(double time, double dt,
                                                               Core::Vector positionBgMolecule,
                                                               Core::Vector velocityBgMolecule,
                                                               Core::Vector positionMolecule, std::vector<Core::Vector> forceMolecules,
                                                               double distance, bool endOfTrajectory) {
    if (time>nextTrajectorySampleTime_) {
        nextTrajectorySampleTime_ += minimalSampleInterval_;
        *trajectoryOutputStream_ << positionBgMolecule.x() << ", "
                << positionBgMolecule.y() << ", "
                << positionBgMolecule.z() << ", "
                << distance << ", "
                << time << ", "
                << velocityBgMolecule.x() << ", "
                << velocityBgMolecule.y() << ", "
                << velocityBgMolecule.z() << ", "
                << forceMolecules[1].x() << ", "
                << forceMolecules[1].y() << ", "
                << forceMolecules[1].z() << ", "
                << dt << ", "
                << positionMolecule.x() << ", "
                << positionMolecule.y() << ", "
                << positionMolecule.z() <<
                   std::endl;
    }
    if(endOfTrajectory == true){
        writeTrajectoryDelimiter();
    }
}

void CollisionModel::MDTrajectoryWriter::writeTrajectoryDelimiter() {
    *trajectoryOutputStream_ << "###" << std::endl;
}


