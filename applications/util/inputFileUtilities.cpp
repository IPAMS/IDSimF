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

#include "inputFileUtilities.hpp"

std::vector<std::unique_ptr<ParticleSimulation::SimionPotentialArray> >
AppUtils::readPotentialArrayFiles(std::vector<std::string> potentialArrayFilenames,
        std::string basePathStr, double scale, bool fastAdjustPA) {

    double potentialScale = fastAdjustPA ? 1/10000.0 : 1.0;

    std::filesystem::path basePath(basePathStr);
    std::vector<std::unique_ptr<ParticleSimulation::SimionPotentialArray>> potentialArrays;
    std::vector<std::string> potentialArraysNames = potentialArrayFilenames;
    for(const auto &paName: potentialArraysNames){
        std::filesystem::path paPath = basePath / paName;
        std::unique_ptr<ParticleSimulation::SimionPotentialArray> pa_pt =
                std::make_unique<ParticleSimulation::SimionPotentialArray>(paPath, scale, potentialScale);
        potentialArrays.push_back(std::move(pa_pt));
    }

    return potentialArrays;
}




