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

#include "catch.hpp"
#include "Core_randomGenerators.hpp"
#include "CollisionModel_MDInteractionsTrajectorySampler.hpp"
#include "CollisionModel_MDForceField_LJ12_6.hpp"
#include "FileIO_CSVReader.hpp"
#include "FileIO_MolecularStructureReader.hpp"
#include "test_util.hpp"

#include <iostream>


TEST_CASE("Basic test of MD trajectory sampler", "[CollisionModels][MDInteractionsModel]") {
    FileIO::MolecularStructureReader reader = FileIO::MolecularStructureReader();
    std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection = reader.readMolecularStructure("test_molecularstructure_reader.json");
    Core::Particle ion;
    ion.setMolecularStructure(molecularStructureCollection.at("Ar+"));
    ion.setVelocity(Core::Vector(600.0, 50.0, 0.0));
    Core::Vector ionRotation({0,0,0});
    CollisionModel::MDForceField_LJ12_6 forceField(0.205E-30);
    auto forceFieldPtr = std::make_unique<CollisionModel::MDForceField_LJ12_6>(forceField);

    Core::Vector collisionParticlePosition({-50e-10, 1e-10, 0});
    Core::Vector collisionParticleVelocity({1000,0,0});
    Core::Vector collisionParticleRotation;

    CollisionModel::MDInteractionsTrajectorySampler mdSim(
        std::move(forceFieldPtr),
        false,
        molecularStructureCollection,
        nullptr
        );

    mdSim.setLegacyTrajectoryWriter("MD_collisions_trajectory_sampler_test.txt", 10, 0.0, 0);

    //Calculate two trajectories with the same sampler:
    collisionParticleRotation = {0,0,0};
    mdSim.calculateTrajectory(ion, ionRotation,
        "N2", collisionParticlePosition, collisionParticleVelocity, collisionParticleRotation,
        1e-11, 1e-16, 200, true);

    collisionParticleRotation = {0,0,M_PI/2.0};
    mdSim.calculateTrajectory(ion, ionRotation,
        "N2", collisionParticlePosition, collisionParticleVelocity, collisionParticleRotation,
        1e-11, 1e-16, 200, true);

    mdSim.calculateTrajectory(ion, ionRotation,
        "Ar", collisionParticlePosition, collisionParticleVelocity, collisionParticleRotation,
    1e-11, 1e-16, 300, false);

    ion.setMolecularStructure(molecularStructureCollection.at("Acetone"));
    ionRotation = {90,0,0};
    mdSim.calculateTrajectory(ion, ionRotation,
    "Ar", collisionParticlePosition, collisionParticleVelocity, collisionParticleRotation,
1e-11, 1e-16, 300, false);

    FileIO::CSVReader csvReader;
    std::vector<std::vector<std::string>> readBack_result = csvReader.readCSVFile("MD_collisions_trajectory_sampler_test.txt", ',');
    CHECK(readBack_result.size() == 762);
    std::vector<double> times = csvReader.extractDouble(readBack_result, 4);
    std::vector<double> bgMolecule_x = csvReader.extractDouble(readBack_result, 0);

    // Trajectories should begin at the right indices
    CHECK(Approx(times.at(0))==1e-16);
    CHECK(Approx(times.at(200))==1e-16);
    CHECK(Approx(times.at(400))==1e-16);

    // start positions should be equal, but rotation of background molecule should change things
    CHECK(Approx(bgMolecule_x[0])==bgMolecule_x[200]);
    CHECK(Approx(bgMolecule_x[100]) != bgMolecule_x[300]);
}