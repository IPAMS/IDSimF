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
 CollisionModel_util.hpp

 Collection of utility methods / functions for backgrond gas interaction modeling

 ****************************/

#ifndef IONSIMULATION_CPP_COLLISIONMODEL_UTIL_HPP
#define IONSIMULATION_CPP_COLLISIONMODEL_UTIL_HPP

#include "BTree_particle.hpp"
#include "RS_AbstractReaction.hpp"
#include <cmath>
#include <functional>


namespace CollisionModel::util {
    const double M_AIR    = 28.94515; ///<Effective mass of air (amu)
    const double D_AIR    = 0.366;    ///<Effective diameter of air (nm)


    double getAirToGas(double massIon_amu, double diameterIon_nm,
                       double collisionGasMass_amu, double collisionGasDiameter_nm);

    double estimateCollisionDiameterFromMass(double massIon_amu);

    double estimateMobility(double massIon_amu, double diameterIon_nm,
                            double collisionGasMass_amu, double collisionGasDiameter_nm);

    std::function<void(RS::CollisionConditions,BTree::Particle&)> getCollisionCountFunction(int* countVal);
}


#endif //IONSIMULATION_CPP_COLLISIONMODEL_UTIL_HPP
