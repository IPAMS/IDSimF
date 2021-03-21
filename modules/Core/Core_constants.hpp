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
 Core_constants.hpp

 collection of general mathematical / physical / chemical constants

 ****************************/

#ifndef BTree_constants_h
#define BTree_constants_h

#define _USE_MATH_DEFINES

#include <cmath>

namespace Core {
    //typedef std::mt19937 MyRNG;  // the Mersenne Twister with a popular choice of parameters
    
    const double ELEMENTARY_CHARGE = 1.60217e-19; ///< Elementary charge (Coulomb)
    const double EPSILON_0 = 8.854e-12;           ///< Vacuum permittivity (Farad / m)
    const double ELECTRIC_CONSTANT = (4.0*M_PI*EPSILON_0); ///< Electrical constant
    const double AMU_TO_KG = 1.66048e-27;            ///<(kg/amu) conversion factor
    
    const double K_BOLTZMANN = 1.3806505e-23;       ///< Boltzmann constant (J/K)
    const double RGas = 8.3145;                     ///< Ideal gas constant (J/(mol*K))
    const double JOULE_TO_EV  = 6.2415095e+18;        ///< (eV/J) conversion factor
    const double N_AVOGADRO = 6.02214199e23;        ///< Avogadro's number
    const double MOL_VOLUME = 22.413996e-3;         //Volume (m^3) of one mol
                                                    // of ideal gas at 0 C, 101.325 kPa
}
#endif /* BTree_constants_h */
