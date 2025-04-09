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

#include "CollisionModel_MathFunctions.hpp"
#include "Core_randomGenerators.hpp"


/**
 * Generates an uniformly distributed random sample on a sphere with radius 'r'
 */
Core::Vector CollisionModel::sphereRand(double r) {
    // Algorithm references:
    // 1. Marsaglia, G.: Choosing a Point from the Surface of a Sphere. Ann. Math. Statist. 43, 645–646 (1972).
    // https://doi.org/10.1214/aoms/1177692644
    // 2. Knop, R.E.: Algorithm 381: random vectors uniform in solid angle. Commun. ACM. 13, 326 (1970).
    // https://doi.org/10.1145/362349.362377

    double xp;
    double yp;
    double S;

    Core::RandomSource* rndSource = Core::globalRandomGeneratorPool->getThreadRandomSource();
    // Find a sample on the unit disc by rejection sampling:
    do
    {
        xp = 2.0*rndSource->uniformRealRndValue() - 1.0;
        yp = 2.0*rndSource->uniformRealRndValue() - 1.0;
        S = xp*xp + yp*yp;
    }
    while( S > 1.0);

    double z = (2.0*S - 1.0)*r; // The azimuthal position z is uniformly in [-1,1], S is uniformly distributed in [0,1]
    double f = 2.0*r*sqrt(1.0-S); //rescaling factor for x,y
    double x = xp*f;
    double y = yp*f;

    return {x,y,z};
}


double CollisionModel::calcSign(const double value){
    if(value > 0){
        return 1.;
    }else if(value < 0){
        return -1.;
    }else{
        return 0;
    }
}

void CollisionModel::rotate(const Core::Vector &angles, Core::Vector& position){

    double tmp_x = angles.x();
    double tmp_y = angles.y();
    double tmp_z = angles.z();

    double new_rel_x = cos(tmp_y) * cos(tmp_z) * position.x()
                        + (sin(tmp_x) * sin(tmp_y) * cos(tmp_z) + cos(tmp_x) * sin(tmp_z)) * position.y()
                        + (sin(tmp_x) * sin(tmp_z) - cos(tmp_x) * sin(tmp_y) * cos(tmp_z)) * position.z();
    double new_rel_y = - cos(tmp_y) * sin(tmp_z) * position.x()
                        + (cos(tmp_x) * cos(tmp_z) - sin(tmp_x) * sin(tmp_y) * sin(tmp_z)) * position.y()
                        + (cos(tmp_x) * sin(tmp_y) * sin(tmp_z) + sin(tmp_x)*cos(tmp_z)) * position.z();
    double new_rel_z = sin(tmp_y) * position.x()
                        - sin(tmp_x) * cos(tmp_y) * position.y()
                        + cos(tmp_x) * cos(tmp_y) * position.z();

    position.x(new_rel_x);
    position.y(new_rel_y);
    position.z(new_rel_z);
}
