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
 BT-parallelTests.cpp

 Simple Tests of different vector implementations

 ****************************/


#include <iostream>
#include <omp.h>
#include <vector>
#include <ctime>
#include <cmath>
#include "Core_matrix3.hpp"

template<typename T>
std::vector<T> testMatrix3(std::size_t nElements, int nSteps, std::string message){
    clock_t begin = std::clock();
    std::time_t wallBegin = std::time(nullptr);

    std::vector<T> matrices;
    for (std::size_t i = 0; i < nElements; i++) {
        matrices.push_back(T({
            1.0, 2.0, 3.0,
            -1.0, 0.0, -1.0,
            3.0, 2.0, 1.0
            })
            * (0.4999995 + (0.0000001/nElements) * i));
    }

    std::size_t nOperations = 0;
    for (int step = 0; step < nSteps; step++) {
        //#pragma omp parallel for default(none) shared(myVector) private(step, nElements)
        for (std::size_t i = 0; i < nElements - 1; i++) {
            matrices[i] = (matrices[i] * matrices[i+1]);
            nOperations ++;
        }
    }
    std::cout << matrices[0]<< " \n\n"<<matrices[nElements-2]<<std::endl;
    clock_t end = std::clock();
    std::time_t wallEnd = std::time(nullptr);

    double cpu_secs = double(end - begin) / CLOCKS_PER_SEC;
    double wall_secs = double(wallEnd - wallBegin);

    std::cout << message << "\n" <<"operations:"<< nOperations<<", elapsed wall time:"<< wall_secs<<" seconds"<<std::endl;
    std::cout << "cpu time:"<< cpu_secs<<std::endl;
    return(matrices);
}

int main(int argc, const char * argv[]) {
    unsigned int nElements = 40000000;
    int nSteps = 25;

    if (argc <2){
        std::cout << "no mode given"<<std::endl;
        return(0);
    }

    std::string mode = argv[1];

    if (mode == "core"){
        testMatrix3<Core::Matrix3> (nElements, nSteps, "Core");
    }
    return 0;
}
