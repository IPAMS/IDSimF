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
#include "Core_vector.hpp"
#include "catch.hpp"
#include "test_util.hpp"

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

        Core::Matrix3 mat3= {3.0, 2.0, 1.0, -1.0, -2.0, -3.0, 4.0, 5.0, 6.0};
        CHECK(mat3(0,0) == Approx(3.0));
        CHECK(mat3(2,1) == Approx(-3.0));
        CHECK(mat3(1,2) == Approx(5.0));
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

        std::array<double, 9> expected = {1.0, 2.0, 15.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0};
        CHECK(mat1.vectorize() == expected);
    }
}

TEST_CASE("Test Matrix3 operators", "[Core][Matrix3]") {

    Core::Matrix3 matA({1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0});
    Core::Matrix3 matB({2.0, 4.0, 6.0, 3.0, 2.0, 1.0, -6.0, -4.0, -2.0});
    Core::Matrix3 matC({1.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 2.0});
    Core::Matrix3 matU({1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0});

    Core::Vector vecA(1.0, 2.0, 3.0);
    Core::Vector vecB(-3.0, -2.0, -1.0);

    SECTION("Matrix addition and subtraction works") {

        Core::Matrix3 expectedAplusB = {3.0, 6.0, 9.0, 7.0, 7.0, 7.0, 1.0, 4.0, 7.0};
        CHECK_THAT(matA + matB, ApproxEqual(expectedAplusB));
        CHECK_THAT(matB + matA, ApproxEqual(expectedAplusB));

        Core::Matrix3 expectedBplusB = {4.0, 8.0, 12.0, 6.0, 4.0, 2.0, -12.0, -8.0, -4.0};
        CHECK_THAT(matB + matB, ApproxEqual(expectedBplusB));

        Core::Matrix3 expectedAminusB = {-1.0, -2.0, -3.0, 1.0, 3.0, 5.0, 13.0, 12.0, 11.0};
        CHECK_THAT(matA - matB, ApproxEqual(expectedAminusB));

        Core::Matrix3 expectedBminusA = {1.0, 2.0, 3.0, -1.0, -3.0, -5.0, -13.0, -12.0, -11.0};
        CHECK_THAT(matB - matA, ApproxEqual(expectedBminusA));
    }

    SECTION("Matrix * scalar multipliction works") {
        Core::Matrix3 expectedAtimes2 = {2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0};
        CHECK_THAT(matA * 2.0, ApproxEqual(expectedAtimes2));

        Core::Matrix3 expectedCtimesNeg1_5 = {-1.5, 0.0, 0.0, 0.0, -7.5, 0.0, 0.0, 0.0, -3.0};
        CHECK_THAT(matC * -1.5, ApproxEqual(expectedCtimesNeg1_5));
    }

    SECTION("Matrix * Vector multiplication works") {
        Core::Vector expectedAtimesVecA(30.0, 36.0, 42.0);
        CHECK_THAT(matA * vecA, ApproxEqual(expectedAtimesVecA));

        Core::Vector expectedBtimesVecB(-6.0, -12.0, -18.0);
        CHECK_THAT(matB * vecB, ApproxEqual(expectedBtimesVecB));
    }

    SECTION("Matrix * Matrix multiplication works") {
        Core::Matrix3 expectedAtimesB = {60, 72, 84, 18, 24, 30, -36, -48, -60};
        CHECK_THAT(matA * matB, ApproxEqual(expectedAtimesB));

        Core::Matrix3 expectedBtimesA = {-10, -4, 2, -13, 2, 17, -16, 8, 32};
        CHECK_THAT(matB * matA, ApproxEqual(expectedBtimesA));
    }

    SECTION("Matrix equality operator works") {
        Core::Matrix3 matEqual = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
        CHECK(matA == matEqual);

        Core::Matrix3 matDiff1 = {1.1, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
        CHECK(!(matA == matDiff1));

        Core::Matrix3 matDiff2 = {1.0, 2.1, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
        CHECK(!(matA == matDiff2));

        Core::Matrix3 matDiff3 = {1.0, 2.0, 3.1, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
        CHECK(!(matA == matDiff3));

        Core::Matrix3 matDiff4 = {1.0, 2.0, 3.0, 4.1, 5.0, 6.0, 7.0, 8.0, 9.0};
        CHECK(!(matA == matDiff4));

        Core::Matrix3 matDiff5 = {1.0, 2.0, 3.0, 4.0, 5.1, 6.0, 7.0, 8.0, 9.0};
        CHECK(!(matA == matDiff5));

        Core::Matrix3 matDiff6 = {1.0, 2.0, 3.0, 4.0, 5.0, 6.1, 7.0, 8.0, 9.0};
        CHECK(!(matA == matDiff6));

        Core::Matrix3 matDiff7 = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.1, 8.0, 9.0};
        CHECK(!(matA == matDiff7));

        Core::Matrix3 matDiff8 = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.1, 9.0};
        CHECK(!(matA == matDiff8));

        Core::Matrix3 matDiff9 = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.1};
        CHECK(!(matA == matDiff9));
    }

    SECTION("Matrix inequality operator works") {
        Core::Matrix3 matEqual = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
        CHECK(!(matA != matEqual));

        Core::Matrix3 matDiff1 = {1.1, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
        CHECK(matA != matDiff1);

        Core::Matrix3 matDiff2 = {1.0, 2.1, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
        CHECK(matA != matDiff2);

        Core::Matrix3 matDiff3 = {1.0, 2.0, 3.1, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
        CHECK(matA != matDiff3);

        Core::Matrix3 matDiff4 = {1.0, 2.0, 3.0, 4.1, 5.0, 6.0, 7.0, 8.0, 9.0};
        CHECK(matA != matDiff4);

        Core::Matrix3 matDiff5 = {1.0, 2.0, 3.0, 4.0, 5.1, 6.0, 7.0, 8.0, 9.0};
        CHECK(matA != matDiff5);

        Core::Matrix3 matDiff6 = {1.0, 2.0, 3.0, 4.0, 5.0, 6.1, 7.0, 8.0, 9.0};
        CHECK(matA != matDiff6);

        Core::Matrix3 matDiff7 = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.1, 8.0, 9.0};
        CHECK(matA != matDiff7);

        Core::Matrix3 matDiff8 = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.1, 9.0};
        CHECK(matA != matDiff8);

        Core::Matrix3 matDiff9 = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.1};
        CHECK(matA != matDiff9);
    }

    SECTION("Matrix transpose works"){
        Core::Matrix3 mat = {0, 3.0, -2.0, -3.0, 0, 1.0, 2.0, -1.0, 0};
        Core::Matrix3 expectedTranspose = {0, -3.0, 2.0, 3.0, 0, -1.0, -2.0, 1.0, 0};
        CHECK_THAT(mat.transpose(), ApproxEqual(expectedTranspose));
    }

    SECTION("Special matrix operations work") {
        Core::Matrix3 expectedStarVecA = {0, 3.0, -2.0, -3.0, 0, 1.0, 2.0, -1.0, 0};
        CHECK_THAT(Core::star(vecA), ApproxEqual(expectedStarVecA));

        Core::Matrix3 expectedStarVecB = {0, -1.0, 2.0, 1.0, 0, -3.0, -2.0, 3.0, 0};
        CHECK_THAT(Core::star(vecB), ApproxEqual(expectedStarVecB));
    }


}
