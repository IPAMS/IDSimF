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
 md_correctness_test.cpp

 Test of md Correcness

 ****************************/

#include "CollisionModel_MDInteractionsTrajectorySampler.hpp"
#include "CollisionModel_MDForceField_LJ12_6.hpp"
#include "FileIO_MolecularStructureReader.hpp"

#include <iomanip>
#include <iostream>

int main(int argc, char** argv) {
    FileIO::MolecularStructureReader reader = FileIO::MolecularStructureReader();
    std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection = reader.readMolecularStructure("test_molecularstructure_reader.json");
    Core::Particle ion;
    ion.setMolecularStructure(molecularStructureCollection.at("Ar+"));
    ion.setVelocity(Core::Vector(600.0, 50.0, 0.0));
    CollisionModel::MDForceField_LJ12_6 forceField(0.205E-30);
    auto forceFieldPtr = std::make_unique<CollisionModel::MDForceField_LJ12_6>(forceField);

    CollisionModel::ParticleInitialConditions ionInitC{
             {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}
    };

    CollisionModel::ParticleInitialConditions collisionParticleInitialConditions{
             {-50e-10, 1e-10, 0}, {1000,0,0}, {0,0,0}, {0,0,0}
    };

    CollisionModel::MDInteractionsTrajectorySampler mdSim(
     std::move(forceFieldPtr),
     false,
     molecularStructureCollection,
     nullptr
     );

    mdSim.calculateTrajectory(ion, "N2",
ionInitC, collisionParticleInitialConditions,
1e-11, 1e-16, 2, false);
}
