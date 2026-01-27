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
 Matrix3.hpp

 Simple three dimensional square matrix with double type

 ****************************/
#ifndef IDSIMF_CORE_MATRIX3_HPP
#define IDSIMF_CORE_MATRIX3_HPP

#include <fstream>
#include <array>
#include <initializer_list>
#include "Core_vector.hpp"

namespace Core{

    class Matrix3 {

    public:
        // Constructors
        Matrix3() = default;
        explicit Matrix3(const double*);
        Matrix3(std::initializer_list<double>);

        // Access operator:
        double& operator()(const std::size_t row, const std::size_t column);

        // Accessors:
        [[nodiscard]] double element(const std::size_t row, const std::size_t column) const;
        [[nodiscard]] std::array<double, 9> vectorize() const;

        //overloaded operators:
        friend Matrix3 operator+(const Matrix3 &lhs, const Matrix3 &rhs);
        friend Matrix3 operator-(const Matrix3 &lhs, const Matrix3 &rhs);

        // matrix multiplication and vector multiplication:
        friend Matrix3 operator*(const Matrix3 &mat, double scalar);
        friend Vector operator*(const Matrix3 &mat, const Vector &vec);
        friend Matrix3 operator*(const Matrix3 &lhs, const Matrix3 &rhs);


        friend bool operator==(Matrix3 const &lhs, Matrix3 const &rhs);
        friend bool operator!=(const Matrix3 &lhs, const Matrix3 &rhs);

    private:
        double elements_[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    };

    //overloaded operators declarations for namespace:
    Matrix3 operator+(const Matrix3 &lhs, const Matrix3 &rhs);
    Matrix3 operator-(const Matrix3 &lhs, const Matrix3 &rhs);
    Matrix3 operator*(const Matrix3 &mat, double scalar);
    Vector operator*(const Matrix3 &mat, const Vector &vec);
    Matrix3 operator*(const Matrix3 &lhs, const Matrix3 &rhs);
    bool operator==(Matrix3 const &lhs, Matrix3 const &rhs);
    bool operator!=(const Matrix3 &lhs, const Matrix3 &rhs);


    // special operators for molecule rotation calculation
    Matrix3 star(Vector vec);
}

std::ostream& operator <<(std::ostream& out, Core::Matrix3 const& matrix);



#endif //IDSIMF_CORE_MATRIX3_HPP
