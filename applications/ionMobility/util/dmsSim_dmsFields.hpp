/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2022 - Physical and Theoretical Chemistry /
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
 dmsSim_dmsSVFields.hpp

 Field Functions defining DMS Dispersion / Separation Fields

 ****************************/


#ifndef IDSIMF_DMSSIM_DMSFIELDS_HPP
#define IDSIMF_DMSSIM_DMSFIELDS_HPP

#include "appUtils_simulationConfiguration.hpp"
#include "PSim_sampledWaveform.hpp"
#include <functional>

enum CVMode {STATIC_CV, AUTO_CV, MODULATED_CV, MODULATED_AUTO_CV};
enum SVMode {BI_SIN, SQUARE, CLIPPED_SIN};

using SVFieldFctType = std::function<double(double fieldAmplitude_VPerM, double time)>;
using CVFieldFctType = std::function<double(double cvAmplitude_VPerM, double time)>;

SVMode parseSVModeConfiguration(AppUtils::simConf_ptr simConf);
SVFieldFctType createSVFieldFunction(SVMode svMode, double fieldWavePeriod);

CVMode parseCVModeConfiguration(AppUtils::simConf_ptr simConf);
CVFieldFctType createCVFieldFunction(CVMode cvMode, double fieldWavePeriod, AppUtils::simConf_ptr simConf);

#endif //IDSIMF_DMSSIM_DMSFIELDS_HPP
