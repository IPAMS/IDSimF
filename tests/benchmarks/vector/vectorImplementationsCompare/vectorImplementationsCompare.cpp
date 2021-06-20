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
#define EIGEN_NO_DEBUG

#include <iostream>
#include <omp.h>
#include <vector>
#include <ctime>
#include <cmath>
#include <Eigen/Dense>
#include "Core_vector.hpp"

void testParallel(unsigned int nElements, int nSteps){
    std::vector<double> myVector(nElements);
    #pragma omp parallel for
    for (unsigned int i=0; i<nElements; i++){
        myVector[i] = 1e-3*i;
    }


    for (int step=0; step<nSteps; step++) {
    #pragma omp parallel for
        for (unsigned int i = 0; i < nElements; i++) {
            myVector[i] = sin(myVector[i]);
        }
    }
    std::cout <<myVector[0]<< " "<<myVector[500]<<std::endl;
}

template<typename T>
std::vector<T> testVector(unsigned int nElements, int nSteps, std::string message){
    clock_t begin = std::clock();
    std::time_t wallBegin = std::time(nullptr);

    std::vector<T> myVector;
    for (unsigned int i = 0; i < nElements; i++) {
        myVector.push_back(T(1e-3 * i, 2e-3 * i, 3e-3 * i));
    }

    for (int step = 0; step < nSteps; step++) {
        //#pragma omp parallel for default(none) shared(myVector) private(step, nElements)
        for (unsigned int i = 1; i < nElements - 1; i++) {
            myVector[i] = (myVector[i] + myVector[i + 1]) / 2.0;
        }
    }
    std::cout << myVector[0]<< " | "<<myVector[10]<<std::endl;
    clock_t end = std::clock();
    std::time_t wallEnd = std::time(nullptr);

    double cpu_secs = double(end - begin) / CLOCKS_PER_SEC;
    double wall_secs = double(wallEnd - wallBegin);


    std::cout << message << "\n" <<"Elapsed wall time:"<< wall_secs<<" seconds"<<std::endl;
    std::cout << "cpu time:"<< cpu_secs<<std::endl;
    return(myVector);
}

int main(int argc, const char * argv[]) {
    unsigned int nElements = 2000000;
    int nSteps = 1000;

    if (argc <2){
        std::cout << "no mode given"<<std::endl;
        return(0);
    }

    std::string mode = argv[1];

    if (mode == "eigen") {
        testVector<Eigen::Vector3d>(nElements, nSteps, "Eigen");
    }

    if (mode == "core"){
        testVector<Core::Vector>   (nElements, nSteps, "Core");
    }

    if (mode == "all"){
        testVector<Eigen::Vector3d>(nElements, nSteps, "Eigen");
        testVector<Core::Vector>   (nElements, nSteps, "Core");
    }

    return 0;
}
