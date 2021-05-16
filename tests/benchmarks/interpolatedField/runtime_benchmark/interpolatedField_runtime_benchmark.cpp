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

 A simple runtime benchmark for interpolated scalar and vector fields

 ****************************/

#include "PSim_interpolatedField.hpp"
#include "Core_vector.hpp"
#include <iostream>
#include <array>
#include <ctime>

struct testResult {
    double scalarSum;
    Core::Vector vectorSum;
};

template<typename fieldT>
testResult performFieldTest(std::string implementationName, int repeats){

    fieldT scalarField = fieldT("test_linear_scalar_field_01.h5");
    fieldT vectorField = fieldT("test_linear_vector_field_01.h5");

    double sum = 0.0;
    Core::Vector vecSum;

    std::clock_t begin = std::clock();
    for (int i=0; i< repeats; i++){
        sum += scalarField.getInterpolatedScalar(1.0+i * 1.0/repeats, 1.0, 1.0, 0);
        vecSum = vecSum + vectorField.getInterpolatedVector(1.0+i * 1.0/repeats, 1.0, 1.0, 0);
    }
    std::clock_t end = std::clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << implementationName<<": sum:" << sum << " vec:" << vecSum << " elapsed:" << elapsed_secs<<" for "
        << repeats <<" samples"<<std::endl;

    return {sum, vecSum};
}

int main(int argc, const char * argv[]) {

    int n = 40000000;
    testResult resultOwn = performFieldTest<ParticleSimulation::InterpolatedField>("Own implementation", n);
    return 0;
}
