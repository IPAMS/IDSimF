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
 test_MDInteractions.cpp

 Testing of molecular collision model with LJ-12-6 and dipole forces

 ****************************/

#include "CollisionModel_MDInteractionsPreconstructed.hpp"
#include "CollisionModel_Molecule.hpp"
#include "CollisionModel_Atom.hpp"
#include "Core_randomGenerators.hpp"
#include "Core_constants.hpp"
#include "catch.hpp"
#include "FileIO_MolecularStructureReader.hpp"
#include <iostream>


TEST_CASE("Basic test MD energy conservation diatomic", "[CollisionModels][MDInteractionsModel]") {

    Core::globalRandomGeneratorPool = std::make_unique<Core::TestRandomGeneratorPool>();

    double diameterN2 = CollisionModel::MDInteractionsModelPreconstructed::DIAMETER_N2;
    FileIO::MolecularStructureReader reader = FileIO::MolecularStructureReader();
    std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection = reader.readMolecularStructure("test_molecularstructure_reader.csv");
    Core::Particle ion;
    ion.setMolecularStructure(molecularStructureCollection.at("Li+"));
    ion.setVelocity(Core::Vector(600.0, 0.0, 0.0));
    CollisionModel::MDInteractionsModelPreconstructed mdSim = CollisionModel::MDInteractionsModelPreconstructed(2000000, 298, 28, 
                                                                                    diameterN2,
                                                                                    1.7E-30, 
                                                                                    "N2", 
                                                                                    1e-9, 
                                                                                    1E-17, 
                                                                                    2, 1, 
                                                                                    30e-10, 
                                                                                    molecularStructureCollection);

    mdSim.setTrajectoryWriter("MD_collisions_preconstructed_trajectories_newImplementation.txt", 30e-10, 0);
    mdSim.updateModelTimestepParameters(1, 0);
    double dt = 2e-11;
    mdSim.modifyVelocity(ion, dt);


}
