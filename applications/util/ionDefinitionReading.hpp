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

#ifndef IDSIMF_IONDEFINITIONREADING_HPP
#define IDSIMF_IONDEFINITIONREADING_HPP

#include "json.h"
#include "BTree_particle.hpp"
#include <vector>
#include <string>
#include <memory>


namespace AppUtils{

    const std::string ION_CLOUD_FILE_KEY = "ion_cloud_init_file";

    bool isIonCloudDefinitionPresent(Json::Value& confRoot);
}

#endif //IDSIMF_IONDEFINITIONREADING_HPP