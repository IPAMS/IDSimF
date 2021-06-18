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
 test_interpolatedField.cpp

 Testing of interpolated field implementations

 ****************************/

#include "PSim_interpolatedField.hpp"
#include "Core_vector.hpp"
#include "test_util.hpp"
#include <iostream>
#include <vector>
#include <array>

#include "catch.hpp"


TEST_CASE("Test self implemented interpolated field", "[ParticleSimulation][InterpolatedField][file readers]") {

    SECTION("Test interpolated scalar field with hdf5 file") {
        ParticleSimulation::InterpolatedField intField("test_linear_scalar_field_01.h5");

        //std::vector<std::vector<double>> grid = intField.getGrid();

        SECTION("The returned spatial grid positions are correct") {
            std::vector<std::vector<double>> grid = intField.getGrid();

            std::vector<double> gridX = {0, 2, 5, 15};
            CHECK(gridX==grid[0]);
            std::vector<double> gridY = {0, 2, 10};
            REQUIRE(gridY==grid[1]);
            std::vector<double> gridZ = {0, 2, 5, 7, 10};
            REQUIRE(gridZ==grid[2]);
        }

        SECTION("The returned bounds are correct") {
            std::array<double, 6> correct_bounds = {0, 15, 0, 10, 0, 10};
            std::array<double, 6> bounds = intField.getBounds();
            REQUIRE(correct_bounds==bounds);
        }

        SECTION("The returned upper bounds are correct") {
            std::array<std::size_t, 3> upperBounds = {1, 2, 4};
            REQUIRE(upperBounds==intField.findLowerBoundIndices(0.5, 5.0, 7.1));

            std::array<std::size_t, 3> upperBounds2 = {0, 3, 1};
            REQUIRE(upperBounds2==intField.findLowerBoundIndices(-0.5, 50.0, 2.0));
        }

        SECTION("Uninterpolated data can be accessed and is correct") {
            REQUIRE(intField.getScalar(1, 1, 1, 0)==6.0);
            REQUIRE(intField.getScalar(1, 3, 2, 0)==9.0);
            REQUIRE(intField.getScalar(3, 0, 1, 0)==17.0);
            REQUIRE(intField.getScalar(3, 0, 4, 0)==25.0);
        }

        SECTION("Interpolated data can be accessed and is correct") {
            REQUIRE(intField.getInterpolatedScalar(1.0, 1.0, 1.0, 0)==3.0);
            REQUIRE(intField.getInterpolatedScalar(11.5, 5.0, 6.2, 0)==22.7);
            REQUIRE(intField.getInterpolatedScalar(5.0, 5.0, 0.1, 0)==Approx(10.1));
            REQUIRE_THROWS_AS(
                    intField.getInterpolatedScalar(-100.0, -100.0, 0.1, 0)==0.0,
                    std::invalid_argument);
        }

        /*
         * fixme: update tests
        SECTION("Scalar field should throw exception if vector is requested") {
            REQUIRE_THROWS(intField.getInterpolatedVector(0.1, 0.1, 0.1, 0));
        }

        SECTION("Interpolated field (eigen): NAN values in source field should produce NAN in result"){
            //ParticleSimulation::VtkFieldReader reader = ParticleSimulation::VtkFieldReader();
            ParticleSimulation::InterpolatedFieldEigen intFieldWithNANs =
                    ParticleSimulation::InterpolatedFieldEigen(
                            reader.readVTKFile("Sciex_Q0_simple_001_export_field_01.vts"));

            REQUIRE(std::isnan(intFieldWithNANs.getInterpolatedVector(0.01,0.003,0.004, 0).x()) == true);
        }*/
    }

    SECTION("Test interpolated vector field with hdf5 file") {

        ParticleSimulation::InterpolatedField intField =
                ParticleSimulation::InterpolatedField("test_linear_vector_field_01.h5");

        SECTION("Correctly uninterpolated values should be returned") {
            REQUIRE(intField.getVector(0,0,0, 0) == Core::Vector(-10, 5, 1));
            REQUIRE(intField.getVector(3,2,1, 0) == Core::Vector( 23, 5, 1));
            REQUIRE(intField.getVector(3,0,0, 0) == Core::Vector( -7, 5, 1));
        }

        SECTION("Correctly interpolated values should be returned") {
            REQUIRE(intField.getInterpolatedVector(0.01, 0.0, 0.01, 0) == Core::Vector(0.02, 5.0, 1.0));
            REQUIRE(intField.getInterpolatedVector(0.01, 0.0, 0.01, 1) == Core::Vector(0.02, 15.0, 11.0));
            REQUIRE(intField.getInterpolatedVector(7.0, 0.0, 0.1, 0) == Core::Vector(7.1, 5.0, 1.0));
            REQUIRE((intField.getInterpolatedVector(2.0, 0.0, 2.1, 0) - Core::Vector(4.1, 5, 1)).magnitude() < 1e-6);
            REQUIRE_THROWS_AS(
                    intField.getInterpolatedVector(-100.0, -100.0, 0.1, 0) == Core::Vector(0.0, 0.0, 0.0),
                    std::invalid_argument);
        }
    }

    SECTION("Test large interpolated vector field with hdf5 file"){
        ParticleSimulation::InterpolatedField largeVectorField("quad_dev_flow_3d.h5");
        REQUIRE(vectorApproxCompare(
                largeVectorField.getVector(296, 38, 31, 0),
                Core::Vector(3.15582, 0, 0)) == vectorsApproxEqual);
    }
}
