/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2021 - Physical and Theoretical Chemistry /
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
 test_logging.cpp

 Tests of logging facilities

 ****************************/

#include "catch.hpp"
#include "appUtils_logging.hpp"
#include "test_util.hpp"
#include <string>

TEST_CASE( "Test logging", "[ApplicationUtils]") {
    auto logger = AppUtils::createLogger("test.log");
    logger->info("First log message");
    logger->info("Second log message");
    logger->flush();

    // test written contents:
    std::string lfContents = readTextFile("test.log");
    CHECK(lfContents.find("[I] Second log message") != std::string::npos);
}