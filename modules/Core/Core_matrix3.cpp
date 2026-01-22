/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2026 - Physical and Theoretical Chemistry /
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
 Matrix3.cpp

 Description

 ****************************/
#include "Core_matrix3.hpp"

#include <cassert>


Core::Matrix3::Matrix3(const double* elements) {
    elements_[0][0] = elements[0];
    elements_[1][0] = elements[1];
    elements_[2][0] = elements[2];

    elements_[0][1] = elements[3];
    elements_[1][1] = elements[4];
    elements_[2][1] = elements[5];

    elements_[0][2] = elements[6];
    elements_[1][2] = elements[7];
    elements_[2][2] = elements[8];
}

Core::Matrix3::Matrix3(std::initializer_list<double> elements) {
    const double*  elemBegin = elements.begin();
    elements_[0][0] = elemBegin[0];
    elements_[1][0] = elemBegin[1];
    elements_[2][0] = elemBegin[2];
    elements_[0][1] = elemBegin[3];
    elements_[1][1] = elemBegin[4];
    elements_[2][1] = elemBegin[5];
    elements_[0][2] = elemBegin[6];
    elements_[1][2] = elemBegin[7];
    elements_[2][2] = elemBegin[8];
}

double& Core::Matrix3::operator()(std::size_t row, std::size_t column){
    assert(row < 3 && column < 3);
    return elements_[row][column];
}

double Core::Matrix3::element(std::size_t row, std::size_t column) const {
    assert(row < 3 && column < 3);
    return elements_[row][column];
}


std::array<double, 9> Core::Matrix3::vectorize() const {
    std::array<double, 9> result;
    std::size_t linearIndex = 0;
    for (std::size_t j = 0; j < 3; ++j) {
        for (std::size_t i = 0; i < 3; ++i) {
            result[linearIndex] = elements_[i][j];
            ++linearIndex;
        }
    }
    return result;
}

std::ostream& operator<< (std::ostream& os, Core::Matrix3 const& mat){
    os << mat.element(0,0) << ' ' << mat.element(0,1)  << ' ' << mat.element(0,2) << "\n"
       << mat.element(1,0) << ' ' << mat.element(1,1)  << ' ' << mat.element(1,2) << "\n"
       << mat.element(2,0) << ' ' << mat.element(2,1)  << ' ' << mat.element(2,2);
    return os;
}