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

#include "ionDefinitionReading.hpp"
#include "PSim_ionCloudReader.hpp"

/**
 * Returns if a ion cloud definition file was specified in a simulation configuration
 *
 * @param confRoot a simulation configuration
 */
bool AppUtils::isIonCloudDefinitionPresent(Json::Value &confRoot) {
    return confRoot.isMember(AppUtils::ION_CLOUD_FILE_KEY);
}


/**
 * Reads the ions defined by an ion cloud definition file specified by a simulation configuration
 * into particles and particle pointer vectors
 *
 * @param particles vector to read the particles to
 * @param particlePtrs vector of particle pointers to read the particle pointers to
 * @param confRoot a simulation configuration
 */
void AppUtils::readIonDefinitionFromIonCloudFile(
        std::vector<std::unique_ptr<BTree::Particle>>& particles,
        std::vector<BTree::Particle*>& particlePtrs,
        Json::Value& confRoot) {

    std::string ionCloudFileName = confRoot.get("ion_cloud_init_file", 0).asString();
    ParticleSimulation::IonCloudReader reader = ParticleSimulation::IonCloudReader();
    particles = reader.readIonCloud(ionCloudFileName);
    //prepare a vector of raw pointers
    for (const auto& part : particles){
        particlePtrs.push_back(part.get());
    }
}