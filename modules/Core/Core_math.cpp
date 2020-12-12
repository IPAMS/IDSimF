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

#include "Core_math.hpp"


/**
 * Convert from degrees to radians
 */
double Core::degToRad(double phi){
    return( M_PI / 180.0 * phi);
}

/**
 * Convert from radians to degrees
 */
double Core::radToDeg(double phi){
    return( phi / (M_PI/180.0) );
}

/**
 * Convert a cartesian vector 'vec' to polar coordinates
 * (coordinate system convention is y upwards)
 * @return A three dimensional vector with the elements {radius, azimuth, elevation}
 */
Core::Vector Core::cartesianToPolar(Core::Vector vec){
    double r = vec.magnitude();
    double azimuth = -atan2(vec.z(), vec.x());
    double elevation = asin(vec.y() / r);

    return {r, azimuth, elevation};
}

/**
 * Rotates a vector 'vec' an 'angle' around the elevation axis
 * (which is the z axis in cartesian coordinates)
 */
Core::Vector Core::elevationRotate(Core::Vector vec, double angle){
    return {
            cos(angle)*vec.x() - sin(angle)*vec.y(),
            sin(angle)*vec.x() + cos(angle)*vec.y(),
            vec.z()
    };
}

/**
 * Rotates a vector 'vec' an 'angle' around the azimuth axis
 * (which is the y axis in cartesian coordinates)
 */
Core::Vector Core::azimuthRotate(Core::Vector vec, double angle){
    //double phi = degToRad(angle);
    return {
            cos(angle)*vec.x() + sin(angle)*vec.z(),
            vec.y(),
            -sin(angle)*vec.x() + cos(angle)*vec.z()
    };
}