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
 particleDefinitionReading.hpp

 Utility functions for reading / parsion particle / ion definitions from config files

 ****************************/

#ifndef IDSIMF_APPUTILS_IONDEFINITIONREADING_HPP
#define IDSIMF_APPUTILS_IONDEFINITIONREADING_HPP

#include "Core_particle.hpp"
#include "PSim_particleStartZone.hpp"
#include "appUtils_simulationConfiguration.hpp"
#include <vector>
#include <string>
#include <memory>


namespace AppUtils{

    const std::string ION_CLOUD_FILE_KEY = "ion_cloud_init_file";

    bool isIonCloudDefinitionPresent(const SimulationConfiguration& simConf);

    void readIonDefinitionFromIonCloudFile(
            std::vector<std::unique_ptr<Core::Particle>>& particles,
            std::vector<Core::Particle*>& particlePtrs,
            const SimulationConfiguration& simConf);

    std::unique_ptr<ParticleSimulation::ParticleStartZone> getStartZoneFromIonDefinition(const SimulationConfiguration& simConf);

    void setIonsKineticEnergy(
            std::vector<std::unique_ptr<Core::Particle>>& particles,
            const SimulationConfiguration& simConf);

    void readRandomIonDefinition(
            std::vector<std::unique_ptr<Core::Particle>>& particles,
            std::vector<Core::Particle*>& particlePtrs,
            const SimulationConfiguration& simConf);

    void readIonDefinition(
            std::vector<std::unique_ptr<Core::Particle>>& particles,
            std::vector<Core::Particle*>& particlePtrs,
            const SimulationConfiguration& simConf);
}

#endif //IDSIMF_APPUTILS_IONDEFINITIONREADING_HPP
