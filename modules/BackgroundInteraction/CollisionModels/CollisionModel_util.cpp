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

#include "CollisionModel_util.hpp"


/**
 * Calculates a scaling constant airToGas between the reduced ion mobilities of an ion in
 * air K0air and in a collision gas K0gas:
 * K0gas = K0air * airToGas.
 *
 * The ion is characterized by its ion mass (amu) and effective ion diameter (m) while the collision gas to
 * calculate airToGas for is defined by the mass (amu) and effective diameter (m) of collision gas mass particles.
 */
double CollisionModel::util::getAirToGas(double massIon_amu, double diameterIon_nm,
                                         double collisionGasMass_amu, double collisionGasDiameter_nm){

    // Reduced mass with respect to air and gas
    double reduced_mass_air = massIon_amu * M_AIR / (massIon_amu + M_AIR);
    double reduced_mass_gas = massIon_amu * collisionGasMass_amu / (massIon_amu + collisionGasMass_amu);

    // Scaling constant to convert from mobilityAir to mobilityGas.
    double airToGas = pow(((diameterIon_nm + D_AIR) / (diameterIon_nm + collisionGasDiameter_nm) ),2.0) *
                        sqrt(reduced_mass_air / reduced_mass_gas);

    return airToGas;
}

/**
 * Get rough estimate of ion diameter (in nm) from ion mass (in amu).
 * (Emperical formula noted in SDS paper.)
 */
double CollisionModel::util::estimateCollisionDiameterFromMass(double massIon_amu) {
    double dIon_nm = 0.120415405 * pow(massIon_amu,1.0/3.0);
    //double dIon_m = dIon_nm * 1.0e-9;
    return dIon_nm;
}


/**
 * Estimates ion mobility (k0) in (10-4 m2 V-1 s-1) from the mass of an ion (in amu), effective diameter
 * of ion (in nm) and the collision gas particle mass (in amu) and effective diamenter (in nm).
 */
double CollisionModel::util::estimateMobility(double massIon_amu, double diameterIon_nm,
                                              double collisionGasMass_amu, double collisionGasDiameter_nm) {
    //Estimate mobility in air (in 10-4 m2 V-1 s-1) from d_ion
    //Emperical formula noted in paper.

    double logdm = log10(diameterIon_nm);
    double mobilityAir =
            1.0e-5 * pow(10,(4.9137 - 1.4491*logdm - 0.2772*pow(logdm,2.0) + 0.0717*pow(logdm,3.0)));

    //Compute mobilityGas by scaling mobilityAir (in 10-4 m2 V-1 s-1)
    double k0 = mobilityAir * getAirToGas(massIon_amu, diameterIon_nm, collisionGasMass_amu, collisionGasDiameter_nm);

    return k0;
}


std::function<void(RS::CollisionConditions,Core::Particle&)> CollisionModel::util::getCollisionCountFunction(int* countVal){
    return [=](RS::CollisionConditions, Core::Particle& /*ion*/)->void{ (*countVal)++; };
}