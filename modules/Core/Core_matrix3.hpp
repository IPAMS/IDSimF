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

namespace Core{

    class Matrix3 {

    public:
        // Constructors
        Matrix3() = default;
        explicit Matrix3(const double*);
        explicit Matrix3(std::initializer_list<double>);

        // Access operator:
        double& operator()(std::size_t row, std::size_t column);

        // Accessors:
        [[nodiscard]] double element(std::size_t row, std::size_t column) const;
        std::array<double, 9> vectorize() const;

        // Setters:
        //void element(std::size_t row, std::size_t column, double newElement);


    private:
        double elements_[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    };

    // matrix operators:


    bool operator==(Matrix3 const &lhs, Matrix3 const &rhs);
    bool operator!=(const Matrix3 &lhs, const Matrix3 &rhs);



}

std::ostream& operator <<(std::ostream& out, Core::Matrix3 const& matrix);



#endif //IDSIMF_CORE_MATRIX3_HPP
