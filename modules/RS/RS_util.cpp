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

#include "RS_util.hpp"
#include "RS_constants.hpp"
#include "Core_constants.hpp"
#include "Core_randomGenerators.hpp"

/**
 * Generates a Maxwell Boltzmann distributed random velocity sample for a given temperature and particle mass
 *
 * @param temperature_K the temperature of the gas
 * @param gasParticleMass_amu the gas particle mass in amu
 * @return a random sampled velocity vector
 */
Core::Vector RS::util::maxwellBoltzmannRandomVelocity(double temperature_K, double gasParticleMass_amu){

    double  vrStdevGas = std::sqrt(
            RS::K_BOLTZMANN * temperature_K / (gasParticleMass_amu * RS::KG_PER_AMU));

    double vxGas = Core::globalRandomGenerator->normalRealRndValue() * vrStdevGas;
    double vyGas = Core::globalRandomGenerator->normalRealRndValue() * vrStdevGas;
    double vzGas = Core::globalRandomGenerator->normalRealRndValue() * vrStdevGas;
    return {vxGas,vyGas,vzGas};
}