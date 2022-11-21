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
    ion.setVelocity(Core::Vector(600.0, 50.0, 0.0));
    CollisionModel::MDInteractionsModel mdSim = CollisionModel::MDInteractionsModel(2000000, 298, 4.003, 
                                                                                    3*diameterHe,
                                                                                    0.205E-30, 
                                                                                    "He", 
                                                                                    1e-10, 
                                                                                    1E-17, 
                                                                                    1, 4, 
                                                                                    35e-10, 
                                                                                    molecularStructureCollection);

    double dt = 2e-11;
    mdSim.setTrajectoryWriter("MD_collisions_microscopic_trajectories_test.txt", 2e-8, 5);
    mdSim.modifyVelocity(ion, dt);


    CHECK(Approx(ion.getVelocity().x()).margin(0.2) ==  542.9474708306);
    CHECK(Approx(ion.getVelocity().y()).margin(0.2) ==  122.734352068);
    CHECK(Approx(ion.getVelocity().z()).margin(0.2) ==  30.6867665587);


    int timestep = 0;
    double time = 0.0;
    for(int i = 0; i < 4; i++) {
        mdSim.updateModelTimestepParameters(timestep, time);
        mdSim.modifyVelocity(ion, 2e-11);
        timestep++;
        time += dt;
    }

    CHECK(Approx(ion.getVelocity().x()).margin(0.2) ==  581.7963212923);
    CHECK(Approx(ion.getVelocity().y()).margin(0.2) ==  94.9018795506);
    CHECK(Approx(ion.getVelocity().z()).margin(0.4) ==  -41.7490856531);

    std::string readBack_early = readTextFile("MD_collisions_microscopic_trajectories_test.txt");
    CHECK(readBack_early == "");

    for(int i = 0; i < 4; i++) {
        mdSim.updateModelTimestepParameters(timestep, time);
        mdSim.modifyVelocity(ion, 2e-11);
        timestep++;
        time += dt;
    }

    std::ifstream fstream("MD_collisions_microscopic_trajectories_test.txt");
    std::string line;
    long i;
    for (i = 0; std::getline(fstream, line); ++i){
        if (i==100){

            //parse line
            std::string delimiter = ",";
            size_t pos = 0;
            std::string token;
            std::vector<double> values;
             while ((pos = line.find(delimiter)) != std::string::npos){
                token = line.substr(0, pos);
                line.erase(0, pos + delimiter.length());
                values.push_back(std::strtod(token.c_str(), nullptr));
            }  
            // for(auto j : values){
            //     std::cout << j << std::endl;
            // }
            std::vector<double> compareValues = {-3.10337e-10, 2.95449e-09, 9.44405e-10, 
                                                3.11725e-09, 3.50606e-13, -342.899, -1135.13, -374.144, 
                                                2.76963e-19, -2.63676e-18, -8.42841e-19};
            std::vector<double> compareMargins = {2e-10, 2e-9, 2e-9, 2e-9, 2e-7};

            CHECK(values.size() == compareValues.size());

            CHECK(Approx(values[0]).margin(compareMargins[0]) ==  compareValues[0]);
            CHECK(Approx(values[1]).margin(compareMargins[1]) ==  compareValues[1]);
            CHECK(Approx(values[2]).margin(compareMargins[2]) ==  compareValues[2]);
            CHECK(Approx(values[3]).margin(compareMargins[3]) ==  compareValues[3]);
            CHECK(Approx(values[4]).margin(compareMargins[4]) ==  compareValues[4]);

            //compare individual values
        }
    }
    CHECK(i > 1266);
}
