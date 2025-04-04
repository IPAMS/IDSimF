/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2025 - Physical and Theoretical Chemistry /
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
 test_MDInteractionsTrajectorySampler.cpp

 Tests of explicit molecular dynamics trajectory sampler

 ****************************/


#include "CollisionModel_MDInteractionsTrajectorySampler.hpp"
#include "CollisionModel_Molecule.hpp"
#include "CollisionModel_Atom.hpp"
#include "CollisionModel_MDForceField_LJ12_6.hpp"
#include "Core_randomGenerators.hpp"
#include "Core_constants.hpp"
#include "catch.hpp"
#include "FileIO_MolecularStructureReader.hpp"
#include "test_util.hpp"

#include <iostream>


TEST_CASE("Basic test of MD trajectory sampler", "[CollisionModels][MDInteractionsModel]") {
    FileIO::MolecularStructureReader reader = FileIO::MolecularStructureReader();
    std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection = reader.readMolecularStructure("test_molecularstructure_reader.json");
    Core::Particle ion;
    ion.setMolecularStructure(molecularStructureCollection.at("Ar+"));
    ion.setVelocity(Core::Vector(600.0, 50.0, 0.0));
    CollisionModel::MDForceField_LJ12_6 forceField(0.205E-30);
    auto forceFieldPtr = std::make_unique<CollisionModel::MDForceField_LJ12_6>(forceField);

    Core::Vector particlePosition({-50e-10, 1e-10, 0});
    Core::Vector particleVelocity({1000,0,0});

    CollisionModel::MDInteractionsTrajectorySampler mdSim(
        CollisionModel::MDInteractionsTrajectorySampler::DIAMETER_HE,
        "N2",1e-11,1e-16,
        4.0,
        std::move(forceFieldPtr),
        molecularStructureCollection,
        particlePosition,
        particleVelocity);

    mdSim.setTrajectoryWriter("MD_collisions_trajectory_sampler_test.txt", 1.0, 0, 0.0);
    mdSim.updateModelTimestepParameters(1, 0);
    mdSim.modifyVelocity(ion);

    std::string readBack_result = readTextFile("MD_collisions_trajectory_sampler_test.txt");
}