/***************************
Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2024 - Physical and Theoretical Chemistry /
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
 test_SpatialFunctions.cpp

 Testing of spatial functions / fields

 ****************************/

#include "CollisionModel_SpatialFieldFunctions.hpp"
#include "PSim_interpolatedField.hpp"
#include "PSim_simionPotentialArray.hpp"
#include "catch.hpp"

TEST_CASE("Test spatial functions", "[CollisionModels][SpatialFunctions]") {
    SECTION("Constant spatial functions should produce constant values") {
        auto constScalarFunction = CollisionModel::getConstantScalarFunction(10.0);
        CHECK(constScalarFunction({1.0, 2.0, 3.0}) == Approx(10.0));
        CHECK(constScalarFunction({-2.0, -4.0, -6.0}) == Approx(10.0));
    }

    SECTION("Variable spatial function from interpolated field should produce correct values") {
        ParticleSimulation::InterpolatedField intField("test_linear_scalar_field_01.h5");
        CHECK(intField.getInterpolatedScalar(1.0, 1.0, 1.0, 0) == Approx(3.0));


        auto fieldFct = CollisionModel::getVariableScalarFunction(intField);
        CHECK(fieldFct({1.0, 1.0, 1.0}) == Approx(3.0));
        CHECK(fieldFct({1.0, 1.0, 1.0}) == Approx(intField.getInterpolatedScalar(1.0, 1.0, 1.0, 0)));
        CHECK(fieldFct({11.5, 5.0, 6.2}) == Approx(22.7));
        CHECK(fieldFct({5.0, 5.0, 0.1}) == Approx(10.1));
        CHECK_THROWS_AS(
                fieldFct({-100.0, -100.0, 0.1}) == Approx(0.0),
                std::invalid_argument);
    }

    SECTION("Variable spatial vector function from interpolated field should produce correct values") {
        ParticleSimulation::InterpolatedField intField =
        ParticleSimulation::InterpolatedField("test_linear_vector_field_01.h5");
        auto fieldFct = CollisionModel::getVariableVectorFunction(intField);
        auto fieldFct_1 = CollisionModel::getVariableVectorFunction(intField, 1);

        CHECK(fieldFct({0.01, 0.0, 0.01}) == Core::Vector(0.02, 5.0, 1.0));
        CHECK(fieldFct_1({0.01, 0.0, 0.01}) == Core::Vector(0.02, 15.0, 11.0));

        CHECK(fieldFct({7.0, 0.0, 0.1}) == Core::Vector(7.1, 5.0, 1.0));
        CHECK((fieldFct({2.0, 0.0, 2.1}) - Core::Vector(4.1, 5, 1)).magnitude() < 1e-6);
        CHECK_THROWS_AS(
                fieldFct({-100.0, -100.0, 0.1}) == Core::Vector(0.0, 0.0, 0.0),
                std::invalid_argument);
    }

    SECTION("Variable spatial function from mirrored 2D SIMION PA should produce correct values") {
        std::string paFilename = "simion_test_planar_2d.pa";
        ParticleSimulation::SimionPotentialArray simPa(paFilename);

        auto fieldFct = CollisionModel::getVariableScalarFunction(simPa);

        CHECK_THROWS_MATCHES(fieldFct({1.0,2.0,0.0}),
                ParticleSimulation::PotentialArrayException,
                Catch::Matchers::Contains("is not in the planar potential array"));

        CHECK_THROWS_MATCHES(fieldFct({0.1,0.04,0.0}),
                ParticleSimulation::PotentialArrayException,
                Catch::Matchers::Contains("is not in the planar potential array"));

        CHECK(fieldFct({0.099,0.04,0.0}) == Approx(19.79632));
        CHECK(fieldFct({0.03, 0.02, 0.0}) == Approx(simPa.getInterpolatedPotential(0.03, 0.02, 0.0)));
    }

    SECTION("Variable spatial function from 3D SIMION PA should produce correct values") {
        std::string paFilename = "simion_test_planar_3d.pa";
        ParticleSimulation::SimionPotentialArray simPa(paFilename);

        auto fieldFct = CollisionModel::getVariableScalarFunction(simPa);


        CHECK_THROWS_MATCHES(fieldFct({1.0,2.0,1.0}),
                ParticleSimulation::PotentialArrayException,
                Catch::Matchers::Contains("not in the planar potential array"));

        CHECK(fieldFct({0.02, 0.015, 0.03}) == Approx(-85.37619));
        CHECK(fieldFct({0.038, 0.018, 0.038}) == Approx(10.0));
        CHECK(fieldFct({0.032, 0.008, 0.018}) == Approx(simPa.getInterpolatedPotential(0.032, 0.008, 0.018)));
    }
}
