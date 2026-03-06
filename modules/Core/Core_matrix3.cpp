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

double& Core::Matrix3::operator()(const std::size_t row, const std::size_t column){
    assert(row < 3 && column < 3);
    return elements_[row][column];
}

Core::Vector Core::Matrix3::column(const std::size_t column){
    return Core::Vector({elements_[0][column], elements_[1][column], elements_[2][column]});
}

double Core::Matrix3::element(const std::size_t row, const std::size_t column) const {
    assert(row < 3 && column < 3);
    return elements_[row][column];
}

/**
 * Vectorizes matrix to a linear vector
 */
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

/**
 * Element wise matrix sum
 */
Core::Matrix3 Core::operator+(const Matrix3& lhs, const Matrix3& rhs) {
    return {
            lhs.elements_[0][0]+rhs.elements_[0][0], lhs.elements_[1][0]+rhs.elements_[1][0], lhs.elements_[2][0]+rhs.elements_[2][0],
            lhs.elements_[0][1]+rhs.elements_[0][1], lhs.elements_[1][1]+rhs.elements_[1][1], lhs.elements_[2][1]+rhs.elements_[2][1],
            lhs.elements_[0][2]+rhs.elements_[0][2], lhs.elements_[1][2]+rhs.elements_[1][2], lhs.elements_[2][2]+rhs.elements_[2][2]
        };
}

/**
 * Element wise matrix difference
 */
Core::Matrix3 Core::operator-(const Matrix3& lhs, const Matrix3& rhs) {
    return {
            lhs.elements_[0][0]-rhs.elements_[0][0], lhs.elements_[1][0]-rhs.elements_[1][0], lhs.elements_[2][0]-rhs.elements_[2][0],
            lhs.elements_[0][1]-rhs.elements_[0][1], lhs.elements_[1][1]-rhs.elements_[1][1], lhs.elements_[2][1]-rhs.elements_[2][1],
            lhs.elements_[0][2]-rhs.elements_[0][2], lhs.elements_[1][2]-rhs.elements_[1][2], lhs.elements_[2][2]-rhs.elements_[2][2]
        };
}

/**
 * Element wise matrix - scalar product
 */
Core::Matrix3 Core::operator*(const Matrix3& mat, const double scalar) {
    return {
            mat.elements_[0][0]*scalar, mat.elements_[1][0]*scalar, mat.elements_[2][0]*scalar,
            mat.elements_[0][1]*scalar, mat.elements_[1][1]*scalar, mat.elements_[2][1]*scalar,
            mat.elements_[0][2]*scalar, mat.elements_[1][2]*scalar, mat.elements_[2][2]*scalar
        };
}

Core::Vector Core::operator*(const Matrix3& mat, const Vector& vec) {
    return {
        mat.elements_[0][0]*vec.x() + mat.elements_[0][1]*vec.y() + mat.elements_[0][2]*vec.z(),
        mat.elements_[1][0]*vec.x() + mat.elements_[1][1]*vec.y() + mat.elements_[1][2]*vec.z(),
        mat.elements_[2][0]*vec.x() + mat.elements_[2][1]*vec.y() + mat.elements_[2][2]*vec.z() };
}

Core::Matrix3 Core::operator*(const Matrix3& lhs, const Matrix3& rhs) {
    return {
        lhs.elements_[0][0]*rhs.elements_[0][0] + lhs.elements_[0][1]*rhs.elements_[1][0] + lhs.elements_[0][2]*rhs.elements_[2][0],
        lhs.elements_[1][0]*rhs.elements_[0][0] + lhs.elements_[1][1]*rhs.elements_[1][0] + lhs.elements_[1][2]*rhs.elements_[2][0],
        lhs.elements_[2][0]*rhs.elements_[0][0] + lhs.elements_[2][1]*rhs.elements_[1][0] + lhs.elements_[2][2]*rhs.elements_[2][0],

        lhs.elements_[0][0]*rhs.elements_[0][1] + lhs.elements_[0][1]*rhs.elements_[1][1] + lhs.elements_[0][2]*rhs.elements_[2][1],
        lhs.elements_[1][0]*rhs.elements_[0][1] + lhs.elements_[1][1]*rhs.elements_[1][1] + lhs.elements_[1][2]*rhs.elements_[2][1],
        lhs.elements_[2][0]*rhs.elements_[0][1] + lhs.elements_[2][1]*rhs.elements_[1][1] + lhs.elements_[2][2]*rhs.elements_[2][1],

        lhs.elements_[0][0]*rhs.elements_[0][2] + lhs.elements_[0][1]*rhs.elements_[1][2] + lhs.elements_[0][2]*rhs.elements_[2][2],
        lhs.elements_[1][0]*rhs.elements_[0][2] + lhs.elements_[1][1]*rhs.elements_[1][2] + lhs.elements_[1][2]*rhs.elements_[2][2],
        lhs.elements_[2][0]*rhs.elements_[0][2] + lhs.elements_[2][1]*rhs.elements_[1][2] + lhs.elements_[2][2]*rhs.elements_[2][2]
    };
}

// Matrix equality means, exact, floating point equality here, thus deactivate
// floating point comparison warning
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
bool Core::operator==(Matrix3 const& lhs, Matrix3 const& rhs) {
    return (
        lhs.elements_[0][0] == rhs.elements_[0][0] && lhs.elements_[1][0] == rhs.elements_[1][0] && lhs.elements_[2][0] == rhs.elements_[2][0] &&
        lhs.elements_[0][1] == rhs.elements_[0][1] && lhs.elements_[1][1] == rhs.elements_[1][1] && lhs.elements_[2][1] == rhs.elements_[2][1] &&
        lhs.elements_[0][2] == rhs.elements_[0][2] && lhs.elements_[1][2] == rhs.elements_[1][2] && lhs.elements_[2][2] == rhs.elements_[2][2]
    );
}
#pragma GCC diagnostic pop

bool Core::operator!=(const Matrix3& lhs, const Matrix3& rhs) {
    return !(lhs == rhs);
}

Core::Matrix3 Core::Matrix3::transpose() {
    return{
        this->elements_[0][0], this->elements_[0][1], this->elements_[0][2],
        this->elements_[1][0], this->elements_[1][1], this->elements_[1][2],
        this->elements_[2][0], this->elements_[2][1], this->elements_[2][2]
    };
}

void Core::Matrix3::setColumn(const std::size_t column, Core::Vector vec){
    this->elements_[0][column] = vec.x();
    this->elements_[1][column] = vec.y();
    this->elements_[2][column] = vec.z();
}

/**
 * Special operation for molecule (rigid body) rotation calculation according to:
 * An Introduction to Physically Based Modeling: Rigid Body Simulation I—Unconstrained Rigid Body Dynamics
 * David Baraff
 *
 * Returns the matrix
 *     0     -v[2]   v[1]
 *  v[2]        0   -v[0]
 * -v[1]      v[0]     0
 *
 * for an input vector v ={v[0], v[1], v[2]}
 */
Core::Matrix3 Core::star(Vector vec) {
    return{
        0, vec.z(), -vec.y(),
        -vec.z(), 0, vec.x(),
        vec.y(), -vec.x(), 0
    };
}

double Core::norm1(Matrix3 mat){
    double a = mat.element(0,0) + mat.element(1,0) + mat.element(2,0);
    double b = mat.element(0,1) + mat.element(1,1) + mat.element(2,1);
    double c = mat.element(0,2) + mat.element(1,2) + mat.element(2,2);

    return std::max(a, std::max(b,c));

}

std::ostream& operator<< (std::ostream& os, Core::Matrix3 const& mat){
    os << mat.element(0,0) << ' ' << mat.element(0,1)  << ' ' << mat.element(0,2) << "\n"
       << mat.element(1,0) << ' ' << mat.element(1,1)  << ' ' << mat.element(1,2) << "\n"
       << mat.element(2,0) << ' ' << mat.element(2,1)  << ' ' << mat.element(2,2);
    return os;
}

// Modified Gram-Schmidt orthogonalization 
// modified version is more stable in finite-precision math
Core::Matrix3 Core::mgs(Matrix3 mat) {
    Core::Matrix3 result = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    for(int i = 0; i < 3; i++){
        Core::Vector columnNorm = mat.column(i) / mat.column(i).magnitude();
        result.setColumn(i, columnNorm);

        for(int j=i+1; j < 3; j++){
            Core::Vector nextCol = mat.column(j);
            nextCol = nextCol - (columnNorm * nextCol) * columnNorm / (columnNorm * columnNorm);
            mat.setColumn(j, nextCol);
        }
    }
    return result; 
}