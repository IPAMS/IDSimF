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

    CHECK(Approx(ion.getVelocity().x()).margin(0.02) ==  479.5777664095);
    CHECK(Approx(ion.getVelocity().y()).margin(0.02) ==  -45.9512880269);
    CHECK(Approx(ion.getVelocity().z()).margin(0.01) ==  151.2568045226);

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
            do{
                token = line.substr(0, pos);
                line.erase(0, pos + delimiter.length());
                values.push_back(std::strtod(token.c_str(), nullptr));
            }  while ((pos = line.find(delimiter)) != std::string::npos);

            std::vector<double> compareValues = {1.97047e-10, 1.58704e-09, 2.41647e-11, 1.59938e-09, 3.30249e-12};
            std::vector<double> compareMargins = {2e-10, 2e-9, 2e-9, 2e-9, 2e-9};

            CHECK(values.size() == compareValues.size());

            CHECK(Approx(values[0]).margin(compareMargins[0]) ==  compareValues[0]);
            CHECK(Approx(values[1]).margin(compareMargins[1]) ==  compareValues[1]);
            CHECK(Approx(values[2]).margin(compareMargins[2]) ==  compareValues[2]);
            CHECK(Approx(values[3]).margin(compareMargins[3]) ==  compareValues[3]);
            CHECK(Approx(values[4]).margin(compareMargins[4]) ==  compareValues[4]);

            //compare individual values
        }
    }
    CHECK(i > 5000);
}


TEST_CASE("Reproduce wrong collision", "[CollisionModels][MDInteractionsModel]") {

    FileIO::MolecularStructureReader reader = FileIO::MolecularStructureReader();
    std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection = reader.readMolecularStructure("test_molecularstructure_reader.csv");
    CollisionModel::Molecule mole = CollisionModel::Molecule(
            Core::Vector(0.0, 0.0, 0.0),
            Core::Vector(0.0, 0.0, 0.0),
            molecularStructureCollection.at("Ar+"));

    // Construct the background gas particle
    CollisionModel::Molecule bgMole = CollisionModel::Molecule(
            Core::Vector(0.0, 0.0, 0.0),
            Core::Vector(0.0, 0.0, 0.0),
            molecularStructureCollection.at("He"));

    double diameterHe = CollisionModel::MDInteractionsModel::DIAMETER_HE;

    double collisionRadiusScaling = 2.0;

    CollisionModel::MDInteractionsModel mdSim = CollisionModel::MDInteractionsModel(
            2000000, 298, 4.003,
            3*diameterHe,
            0.205E-30,
            "He",
            1e-10,
            1E-16,
            collisionRadiusScaling, 1,
            35e-10,
            molecularStructureCollection);

    mdSim.setTrajectoryWriter("MD_collisions_faulty_collision.txt", 2e-8, 0);
    mdSim.updateModelTimestepParameters(1, 0.0);

    bgMole.setComPos({1.6905050149374716e-09, 7.882860479075753e-10, -2.3496378233986436e-09});
    bgMole.setComVel({-1848.0048078989612, 862.0310727273678, 2556.1253608501447});
    //bgMole.setComVel({-1848.0048078989612, -862.0310727273678, 2556.1253608501447});

    std::vector<CollisionModel::Molecule*> moleculesPtr = {&mole, &bgMole};

    std::cout << bgMole.getComVel()<< std::endl;
    double collisionRadius = collisionRadiusScaling*(mole.getDiameter() + diameterHe)/2.0;
    bool trajectorySuccess = mdSim.rk4InternAdaptiveStep(moleculesPtr, 1e-16, 1e-10, 4*collisionRadius);
    std::cout << bgMole.getComVel()<< std::endl;
    CHECK(Approx(bgMole.getComVel().magnitude()) == 3269.8644041702);
    CHECK(trajectorySuccess);

}
