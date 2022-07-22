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
 test_MolecularStructureReader.cpp

 Testing of molecular structure reader implementation

 ****************************/

#include "FileIO_MolecularStructureReader.hpp"
#include "catch.hpp"
#include "test_util.hpp"
#include "CollisionModel_MolecularStructure.hpp"
#include <iostream>


TEST_CASE("Test molecular structure reader", "[ParticleSimulation][MolecularStructureReader][file readers]") {

    SECTION( "Molecular Structure reader: existing file should open without exception") {
        FileIO::MolecularStructureReader reader = FileIO::MolecularStructureReader();
        std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection;
        REQUIRE_NOTHROW(molecularStructureCollection = reader.readMolecularStructure("test_molecularstructure_reader.csv"));

        auto it1 = molecularStructureCollection.find("Ar2");
        CHECK(it1 != molecularStructureCollection.end());
        auto it2 = molecularStructureCollection.find("He");
        CHECK(it2 != molecularStructureCollection.end());
    }
}