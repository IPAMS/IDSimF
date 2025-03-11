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
 Core_math.hpp

 Mathematical utility functions of general purpose

 ****************************/
#ifndef IDSIMF_CORE_MATH_HPP
#define IDSIMF_CORE_MATH_HPP

#include "Core_vector.hpp"
#include "Core_constants.hpp"
#include <functional>
#include <iostream>

namespace Core {

    double degToRad(double phi);
    double radToDeg(double phi);

    template<typename T>
    double goldenSectionSearch(const T& objA, const T& objB, std::function<double (const T&, const T&, double)> func,
                                double right, double left, double tol = 1e-5){
    double invphi = (sqrt(5) - 1) / 2;
    while((left-right) > tol){
        double c = left - (left - right) * invphi;
        double d = right + (left - right) * invphi;
        if(func(objA, objB, c) > func(objA, objB, d)){
            left = d;
        }else{
            right = c;
        }
    }
    return (right+left)/2;
}

    Core::Vector cartesianToPolar(Core::Vector vec);
    Core::Vector elevationRotate(Core::Vector vec, double angle);
    Core::Vector azimuthRotate(Core::Vector vec, double angle);


}


#endif //IDSIMF_CORE_MATH_HPP
