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
 test_matrix3.cpp

 Test of simple 3x3 matrix class

 ****************************/

#include "Core_matrix3.hpp"
#include "catch.hpp"
#include "test_util.hpp"
#include "../../../libs/CLI11/CLI11.hpp"

TEST_CASE("Test Matrix3 construction and basic element access", "[Core][Matrix3]") {
    SECTION("Matrix initialization works"){
        Core::Matrix3 zeroMat;
        CHECK(zeroMat(0,0) == Approx(0.0));
        CHECK(zeroMat(1,2) == Approx(0.0));

        Core::Matrix3 mat1({1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0});
        CHECK(mat1(0,0) == Approx(1.0));
        CHECK(mat1(2,0) == Approx(3.0));
        CHECK(mat1(1,2) == Approx(8.0));

        double elements[] = {1.0, 2.0, 3.0, 0.0, 0.0, 0.0, 5.0, 5.0, 5.0};
        Core::Matrix3 mat2(elements);
        CHECK(mat2(0,0) == Approx(1.0));
        CHECK(mat2(1,2) == Approx(5.0));

    }

    SECTION("Element access and assignment works"){
        Core::Matrix3 mat1({1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0});
        CHECK(mat1(2,0) == Approx(3.0));

        mat1(2,0) = 15.0;
        CHECK(mat1(2,0) == Approx(15.0));

        mat1(2,2) = 10.0;
        CHECK(mat1(2,2) == Approx(10.0));


        Core::Matrix3 mat2({1.0, 2.0, 15.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0});
        CHECK_THAT(mat1, ApproxEqual(mat2));

        Core::Matrix3 mat3({2.0, 0.0, 0.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0});
        CHECK_THAT(mat1, !ApproxEqual(mat3));
    }
}
