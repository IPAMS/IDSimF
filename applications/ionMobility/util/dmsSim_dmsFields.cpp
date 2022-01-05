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
 ****************************/

#include "dmsSim_dmsFields.hpp"
#include <cmath>

SVMode parseSVModeConfiguration(AppUtils::simConf_ptr simConf){
    SVMode svMode;
    std::string svModeStr = simConf->stringParameter("sv_mode");
    if (svModeStr == "bi_sin"){
        svMode = BI_SIN;
    }
    else if (svModeStr == "square"){
        svMode = SQUARE;
    }
    else if (svModeStr == "clipped_sin"){
        svMode = CLIPPED_SIN;
    }
    else{
        throw std::invalid_argument("wrong configuration value: sv_mode");
    }
    return svMode;
}

SVFieldFctType createSVFieldFunction(SVMode svMode, double fieldWavePeriod){

    SVFieldFctType SVFieldFct;

    if (svMode == BI_SIN){
        double field_h = 2.0;
        double field_F = 2.0;
        double field_W = (1/fieldWavePeriod) * 2 * M_PI;
        auto  fieldFctBisinusoidal=
                [field_F, field_W, field_h]
                        (double svAmplitude_VPerM, double time) -> double{

                    //double particleCharge = particle->getCharge();
                    double voltageSVgp = svAmplitude_VPerM * 0.6667; // V/m (1V/m peak to peak is 0.6667V/m ground to peak)
                    double voltageSVt = (field_F * sin(field_W * time)
                            + sin(field_h * field_W * time - 0.5 * M_PI))
                            * voltageSVgp / (field_F + 1);

                    return voltageSVt;
                };
        SVFieldFct = fieldFctBisinusoidal;
    }
    else if (svMode == SQUARE){
        double thirdOfWavePeriod = fieldWavePeriod / 3.0;
        auto  fieldFctSquare=
                [fieldWavePeriod, thirdOfWavePeriod]
                        (double svAmplitude_VPerM, double time) -> double{

                    double timeInPeriod = std::fmod(time, fieldWavePeriod);

                    double voltageSVgp_highField = svAmplitude_VPerM * 0.666667; // V/m (1V/m peak to peak is 0.6667V/m ground to peak)
                    double voltageSVgp_lowField = -voltageSVgp_highField * 0.5; // low field is 1/2 of high field
                    if (timeInPeriod < thirdOfWavePeriod){
                        return voltageSVgp_highField;
                    }
                    else {
                        return voltageSVgp_lowField;
                    }
                };
        SVFieldFct = fieldFctSquare;
    }
    else if (svMode == CLIPPED_SIN){
        double h = 3.0/2.0* M_PI  - 1.0;
        double t_sin = M_PI / (2*h + 1);
        double f_low = -(2.0* t_sin) / (M_PI-2.0*t_sin);

        auto  fieldFctClippedSin=
                [f_low, t_sin, fieldWavePeriod]
                        (double svAmplitude_VPerM, double time) -> double{

                    double normalizedTimeInPeriod = std::fmod(time, fieldWavePeriod) / fieldWavePeriod;
                    if (normalizedTimeInPeriod < t_sin){
                        return  ((M_PI * sin(M_PI* normalizedTimeInPeriod / t_sin) - 2*t_sin) / (M_PI-2*t_sin) * svAmplitude_VPerM);
                    }
                    else {
                        return  (f_low*svAmplitude_VPerM);
                    }
                };
        SVFieldFct = fieldFctClippedSin;
    }

    return SVFieldFct;
}

CVMode parseCVModeConfiguration(AppUtils::simConf_ptr simConf){
    std::string cvModeStr = simConf->stringParameter("cv_mode");
    CVMode cvMode;
    if (cvModeStr == "static"){
        cvMode = STATIC_CV;
    }
    else if (cvModeStr == "auto"){
        cvMode = AUTO_CV;
    }
    else if (cvModeStr == "modulated"){
        cvMode = MODULATED_CV;
    }
    else if (cvModeStr == "modulated_auto"){
        cvMode = MODULATED_AUTO_CV;
    }
    else{
        throw std::invalid_argument("wrong configuration value: cv_mode");
    }

    return cvMode;
}
