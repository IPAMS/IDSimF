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

#include "CollisionModel_MDInteractions.hpp"
#include "CollisionModel_Molecule.hpp"
#include "CollisionModel_Atom.hpp"
#include "Core_randomGenerators.hpp"
#include "Core_constants.hpp"
#include "catch.hpp"
#include "FileIO_MolecularStructureReader.hpp"
#include <iostream>

std::string readTextFile(std::string filename){
    std::ifstream ifs(filename);
    std::string content;
    content.assign( (std::istreambuf_iterator<char>(ifs) ), (std::istreambuf_iterator<char>() ) );

    return content;
}

TEST_CASE("Basic test MD Interactions model", "[CollisionModels][MDInteractionsModel]") {

    Core::globalRandomGeneratorPool = std::make_unique<Core::TestRandomGeneratorPool>();

    double diameterHe = CollisionModel::MDInteractionsModel::DIAMETER_HE;
    FileIO::MolecularStructureReader reader = FileIO::MolecularStructureReader();
    std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection = reader.readMolecularStructure("test_molecularstructure_reader.csv");
    Core::Particle ion;
    ion.setMolecularStructure(molecularStructureCollection.at("Ar+"));
    ion.setVelocity(Core::Vector(600.0, 0.0, 0.0));
    CollisionModel::MDInteractionsModel mdSim = CollisionModel::MDInteractionsModel(2000000, 298, 4.003, 
                                                                                    3*diameterHe,
                                                                                    0.205E-30, 
                                                                                    "He", 
                                                                                    1e-10, 
                                                                                    1E-16, 
                                                                                    1, 4, 
                                                                                    35e-10, 
                                                                                    molecularStructureCollection);

    double dt = 2e-11;
    mdSim.setTrajectoryWriter("MD_collisions_microscopic_trajectories_test.txt", 2e-8, 5);
    mdSim.modifyVelocity(ion, dt);

    CHECK(Approx(ion.getVelocity().x()) ==  616.7784099091);
    CHECK(Approx(ion.getVelocity().y()) ==  9.9825566288);
    CHECK(Approx(ion.getVelocity().z()) ==  -19.0017715144);


    int timestep = 0;
    double time = 0.0;
    for(int i = 0; i < 4; i++) {
        mdSim.updateModelTimestepParameters(timestep, time);
        mdSim.modifyVelocity(ion, 2e-11);
        timestep++;
        time += dt;
    }

    std::string readBack_early = readTextFile("MD_collisions_microscopic_trajectories_test.txt");
    CHECK(readBack_early == "");

    for(int i = 0; i < 4; i++) {
        mdSim.updateModelTimestepParameters(timestep, time);
        mdSim.modifyVelocity(ion, 2e-11);
        timestep++;
        time += dt;
    }

    std::string readBack_late = readTextFile("MD_collisions_microscopic_trajectories_test.txt");
    CHECK(readBack_late.size() == 611638);
}
