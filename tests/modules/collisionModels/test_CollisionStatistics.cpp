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
 test_CollisionStatistic.cpp

 Testing of collision statistics for SDS collision modeling

 ****************************/

#include "catch.hpp"
#include "CollisionModel_CollisionStatistics.hpp"

TEST_CASE( "Test basic collision statistics semantics", "[CollisionModels][CollisionStatistics]") {

    SECTION("CollisionStatistics should trow on construction non existing file") {
        REQUIRE_THROWS(CollisionModel::CollisionStatistics("FileNonExisting.dat"));
    }

    SECTION("Test default CollisionStatistics "){
        CollisionModel::CollisionStatistics cs;

        SECTION("CollisionStatistics should initialized with default statistics if no file is given") {
            REQUIRE(cs.getNCollisions() == 100000);
            REQUIRE(cs.getICDFs()[0][6] == Approx(111.8694458));
        }

        SECTION("CollisionStatistics statistics index lookup should work correctly", "[Collision::CollisionStatistics]") {
            REQUIRE(cs.findUpperDistIndex(log10(1.5)) == 0);
            REQUIRE(cs.findUpperDistIndex(log10(9.999)) == 0);
            REQUIRE(cs.findUpperDistIndex(log10(10)) == 1);
            REQUIRE(cs.findUpperDistIndex(log10(20)) == 1);
            REQUIRE(cs.findUpperDistIndex(log10(100.1)) == 2);
            REQUIRE(cs.findUpperDistIndex(log10(5000)) == 3);
            REQUIRE(cs.findUpperDistIndex(log10(100000)) == 4);
            REQUIRE(cs.findUpperDistIndex(log10(10000000)) == 4);
        }
    }

    SECTION("CollisionStatistics should be constructable with existing statistics file correctly"){
        CollisionModel::CollisionStatistics cs("cs_icdf_2020_02_27_001_test.dat");
        REQUIRE(cs.getNDist() == 5);
        REQUIRE(cs.getNDistPoints() == 1002);
        REQUIRE(cs.getNCollisions() == 10000);

        std::vector<std::vector<double>> icdfs = cs.getICDFs();
        REQUIRE(icdfs.size() == 5);
        REQUIRE(Approx(icdfs[0][0]) == 22.104078);
        REQUIRE(Approx(icdfs[2][10]) == 299.96972);

        std::vector<double> testMassRatios = {1,10,100,1000,10000};
        REQUIRE(cs.getMassRatios() == testMassRatios);
    }
}