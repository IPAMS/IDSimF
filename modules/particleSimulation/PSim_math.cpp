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

#include "PSim_math.hpp"

/**
 * Computes an linearly spaced vector of double values (similar to the corresponding matlab and numpy methods)
 * (The returned vector contains the mathematical interval [lower,upper], so the boundary values are
 * both part of the returned vector)
 *
 * @param lower a lower boundary of the computed values
 * @param upper a upper boundary (which
 * @param n the number of calculated values
 * @returns a std::vector with linearly spaced values between [lower,upper]
 */
std::vector<double> ParticleSimulation::linspace(double lower, double upper, int n) {
    std::vector<double> result = std::vector<double>();
    double dif = (upper-lower) / (n-1);
    for (int i=0; i<n; ++i){
        result.push_back(lower+i*dif);
    }
    return(result);
}

/**
 * Fill a vector with a given value
 *
 * @param value the value to fill in the vector
 * @param n length of the vector to create
 * @return a std::vector of length n filled with the value
 */
std::vector<double> ParticleSimulation::fillVector(double value, int n) {
    std::vector<double> result = std::vector<double>();
    for (int i=0; i<n; ++i){
        result.push_back(value);
    }
    return(result);
}