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
 RS_constants.hpp

 Some physical constants

 ****************************/

#ifndef RS_constants_hpp
#define RS_constants_hpp

namespace RS {
    //typedef std::mt19937 MyRNG;  // the Mersenne Twister with a popular choice of parameters

    const double T_standard = 298.15; ///< the standard temperature (K)
    const double kBoltzmann = 1.3806505e-23;       //Boltzmann constant (J/K)
    const double RGas = 8.3145;                    //Ideal gas constant (J/(mol*K))
    const double kgPerAmu = 1.66053873e-27;        //mass per atomic mass unit
}


#endif //RS_constants_hpp
