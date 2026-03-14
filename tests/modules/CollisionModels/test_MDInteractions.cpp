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
#include "CollisionModel_MDForceField_LJ12_6.hpp"
#include "Core_randomGenerators.hpp"
#include "catch.hpp"
#include "FileIO_MolecularStructureReader.hpp"
#include "test_hdf5_util.hpp"
#include <iostream>

std::string readTextFile(std::string filename){
    std::ifstream ifs(filename);
    std::string content;
    content.assign( (std::istreambuf_iterator<char>(ifs) ), (std::istreambuf_iterator<char>() ) );
    return content;
}

TEST_CASE("Basic test MD Interactions model", "[CollisionModels][MDInteractionsModel]") {

    Core::globalRandomGeneratorPool = std::make_unique<Core::XoshiroTestRandomGeneratorPool>();

    double diameterHe = CollisionModel::MDInteractionsModel::DIAMETER_HE;
    FileIO::MolecularStructureReader reader = FileIO::MolecularStructureReader();
    std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection = reader.readMolecularStructure("test_molecularstructure_reader.json");
    Core::Particle ion;
    ion.setMolecularStructure(molecularStructureCollection.at("Ar+"));
    ion.setVelocity(Core::Vector(600.0, 50.0, 0.0));
    CollisionModel::MDForceField_LJ12_6 forceField(0.205E-30);
    auto forceFieldPtr = std::make_unique<CollisionModel::MDForceField_LJ12_6>(forceField);
    CollisionModel::MDInteractionsModel mdSim = CollisionModel::MDInteractionsModel(2000000, 298,
                                                                                    4.003,
                                                                                    diameterHe,
                                                                                    "He",
                                                                                    1e-10, 
                                                                                    1E-17,
                                                                                    2, 1,
                                                                                    35e-10,
                                                                                    std::move(forceFieldPtr),
                                                                                    molecularStructureCollection);

    double dt = 2e-11;

    SECTION("Test basic MD Model with HDF5 trajectory writer") {
        mdSim.setHDF5TrajectoryWriter("MD_collisions_resolved_trajectories_test.h5", 35e-10, 0);
        mdSim.modifyVelocity(ion, dt);
        CHECK(Approx(ion.getVelocity().x()).margin(0.2) ==  449.2092547232);
        CHECK(Approx(ion.getVelocity().y()).margin(0.2) ==  -36.8772475434);
        CHECK(Approx(ion.getVelocity().z()).margin(0.2) ==  45.5651248115);

        unsigned int timestep = 0;
        double time = 0.0;
        for(int i = 0; i < 4; i++) {
            mdSim.updateModelTimestepParameters(timestep, time);
            mdSim.modifyVelocity(ion, 2e-11);
        }

        CHECK(Approx(ion.getVelocity().x()).margin(0.8) ==  252.9988351158); //252.9988351158);
        CHECK(Approx(ion.getVelocity().y()).margin(0.8) ==  -170.992193862); //-170.992193862);
        CHECK(Approx(ion.getVelocity().z()).margin(0.2) ==  -267.150091929); //-267.150091929);
    }


    SECTION("Test basic MD Model with legacy trajectory writer") {
        mdSim.setLegacyTrajectoryWriter("MD_collisions_microscopic_trajectories_test.txt", 35e-10, 0);
        mdSim.modifyVelocity(ion, dt);
        std::cout << ion.getVelocity() << std::endl;


        CHECK(Approx(ion.getVelocity().x()).margin(0.2) ==  449.2092547232);
        CHECK(Approx(ion.getVelocity().y()).margin(0.2) ==  -36.8772475434);
        CHECK(Approx(ion.getVelocity().z()).margin(0.2) ==  45.5651248115);


        unsigned int timestep = 0;
        double time = 0.0;
        for(int i = 0; i < 4; i++) {
            mdSim.updateModelTimestepParameters(timestep, time);
            mdSim.modifyVelocity(ion, 2e-11);
        }

        CHECK(Approx(ion.getVelocity().x()).margin(0.8) ==  252.9988351158); //252.9988351158);
        CHECK(Approx(ion.getVelocity().y()).margin(0.8) ==  -170.992193862); //-170.992193862);
        CHECK(Approx(ion.getVelocity().z()).margin(0.2) ==  -267.150091929); //-267.150091929);

        for(int i = 0; i < 4; i++) {
            mdSim.updateModelTimestepParameters(timestep, time);
            mdSim.modifyVelocity(ion, 2e-11);
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
                //parse last value
                if (line.size() > 0) {
                    values.push_back(std::strtod(line.c_str(), nullptr));
                }
                // for(auto j : values){
                //     std::cout << j << std::endl;
                // }
                std::vector<double> compareValues =
                    {1.34745e-10, -2.27634e-10, 3.09759e-11, 2.6301e-10, 2.46681e-12, -1379.11, -28.5662, -112.195, 1.09982e-10, -1.87699e-10, 2.54482e-11, 1.16555e-15, 3.47231e-12, -2.22873e-12, 4.82311e-13};
                    //{1.3633e-10, -2.2762e-10, 3.11065e-11, 2.64082e-10, 2.46567e-12, -1397.99, 3.66005, -116.564, 1.03247e-10, -1.73467e-10, 2.3652e-11, 1.14127e-15, 3.46571e-12, -2.24046e-12, 4.82781e-13};
                /*std::vector<double> compareValues = {0.10036e-10, 0.62538e-09, -2.34314e-10, 2.78511e-10,
                                                    6.47741e-13, -342.483, -1132.27, -375.806, 8.4005e-17,
                                                    -5.37867e-16, -1.70928e-16};*/
                std::vector<double> compareMargins = {2e-10, 2e-9, 2e-9, 2e-9, 2e-7};

                CHECK(values.size() == compareValues.size());

                CHECK(Approx(values[0]).margin(compareMargins[0]) ==  compareValues[0]);
                // CHECK(Approx(values[1]).margin(compareMargins[1]) ==  compareValues[1]);
                CHECK(Approx(values[2]).margin(compareMargins[2]) ==  compareValues[2]);
                CHECK(Approx(values[3]).margin(compareMargins[3]) ==  compareValues[3]);
                CHECK(Approx(values[4]).margin(compareMargins[4]) ==  compareValues[4]);

                //compare individual values
            }
        }
        CHECK(i > 920);
    }
}

bool checkDims(std::string filename, std::string dataSetName, hsize_t xDim, hsize_t yDim) {
    auto tra1DS = openDataSet(filename.c_str(), dataSetName.c_str());
    hsize_t dims[2];
    hsize_t maxDims[2];
    H5Sget_simple_extent_dims(tra1DS.getSpace().getId(), dims, maxDims);

    if (dims[0] == xDim && dims[1] == yDim && maxDims[0] == H5S_UNLIMITED && maxDims[1] == yDim) {
        return true;
    }
    UNSCOPED_INFO("check dims failed! dims[0]: "<<dims[0]<<" dims[1]: "<<dims[1]<<" xDim: "<<xDim <<" yDim: "<<yDim);
    return false;
}

TEST_CASE("Test MD Model with Multi-Atom Molecules", "[CollisionModels][MDInteractionsModel]") {
    Core::globalRandomGeneratorPool = std::make_unique<Core::XoshiroTestRandomGeneratorPool>();

    FileIO::MolecularStructureReader reader = FileIO::MolecularStructureReader();
    std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection =
        reader.readMolecularStructure("test_molecularstructure_reader.json");
    Core::Particle ion;
    ion.setMolecularStructure(molecularStructureCollection.at("O2+"));
    ion.setVelocity(Core::Vector(600.0, 50.0, 0.0));
    CollisionModel::MDForceField_LJ12_6 forceField(0.205E-30);
    auto forceFieldPtr = std::make_unique<CollisionModel::MDForceField_LJ12_6>(forceField);

    double dt = 2e-11;

    /*SECTION("Test with He as collision gas") {
        CollisionModel::MDInteractionsModel mdSim = CollisionModel::MDInteractionsModel(
            2000000, 298,
            4.003, CollisionModel::MDInteractionsModel::DIAMETER_HE, "He",
            1e-10, 1E-17, 2, 1, 35e-10,
            std::move(forceFieldPtr), molecularStructureCollection);

        std::string h5Filename = "MD_collisions_multiatom_trajectories_He.h5";
        mdSim.setHDF5TrajectoryWriter(h5Filename, 35e-10, 0);
        mdSim.modifyVelocity(ion, dt);
        unsigned int timestep = 0;
        double time = 0.0;
        for(int i = 0; i < 4; i++) {
            mdSim.updateModelTimestepParameters(timestep, time);
            mdSim.modifyVelocity(ion, 2e-11);
        }

        // check results written into trajectory:
        CHECK(checkDims(h5Filename, "MD_trajectories/trajectory1", 197, 11));
        CHECK(checkDims(h5Filename, "MD_trajectories/trajectory2", 193, 11));
        H5::Group group = openGroup(h5Filename, "MD_trajectories");
        CHECK(group.getNumObjs() == 3);

        //check Attributes of first trajectory:
        H5::DataSet tra1DS = openDataSet(h5Filename, "MD_trajectories/trajectory1");
        auto colNamesAttribute = readAttribute<H5::DataSet, std::string>(tra1DS, "column_names");
        std::vector<std::string> expectedColNames = {"time", "dt",
            "mol_a0_pos_x", "mol_a0_pos_y", "mol_a0_pos_z", "mol_a1_pos_x", "mol_a1_pos_y", "mol_a1_pos_z",
            "bg_a0_pos_x", "bg_a0_pos_y", "bg_a0_pos_z"};
        CHECK_THAT(colNamesAttribute, Catch::Matchers::Equals(expectedColNames));

        auto nAtomsAttribute = readAttribute<H5::DataSet, std::size_t>(tra1DS, "number_of_atoms");
        std::vector<std::size_t> expectedAtomNumbers = {2, 1};
        CHECK_THAT(nAtomsAttribute, Catch::Matchers::Equals(expectedAtomNumbers));
    }*/

    SECTION("Test with N2 as collision gas") {
        CollisionModel::MDInteractionsModel mdSim = CollisionModel::MDInteractionsModel(
            2000000, 298,
            28.0, CollisionModel::MDInteractionsModel::DIAMETER_N2, "N2",
            1e-10, 1E-17, 2, 1, 35e-10,
            std::move(forceFieldPtr), molecularStructureCollection);

        std::string h5Filename = "MD_collisions_multiatom_trajectories_N2.h5";
        mdSim.setHDF5TrajectoryWriter(h5Filename, 35e-10, 0);
        mdSim.modifyVelocity(ion, dt);
        unsigned int timestep = 0;
        double time = 0.0;
        for(int i = 0; i < 4; i++) {
            mdSim.updateModelTimestepParameters(timestep, time);
            mdSim.modifyVelocity(ion, 2e-11);
        }

        // check results written into trajectory:
        CHECK(checkDims(h5Filename, "MD_trajectories/trajectory1", 135, 17));
        CHECK(checkDims(h5Filename, "MD_trajectories/trajectory2", 183, 17));
        H5::Group group = openGroup(h5Filename, "MD_trajectories");
        CHECK(group.getNumObjs() == 4);

        //check Attributes of first trajectory:
        H5::DataSet tra1DS = openDataSet(h5Filename, "MD_trajectories/trajectory1");
        auto colNamesAttribute = readAttribute<H5::DataSet, std::string>(tra1DS, "column_names");
        std::vector<std::string> expectedColNames = {"time", "dt",
            "mol_a0_pos_x", "mol_a0_pos_y", "mol_a0_pos_z", "mol_a1_pos_x", "mol_a1_pos_y", "mol_a1_pos_z",
            "bg_a0_pos_x", "bg_a0_pos_y", "bg_a0_pos_z", "bg_a1_pos_x", "bg_a1_pos_y", "bg_a1_pos_z", "bg_a2_pos_x", "bg_a2_pos_y", "bg_a2_pos_z"};
        CHECK_THAT(colNamesAttribute, Catch::Matchers::Equals(expectedColNames));

        auto atomMasses = readAttribute<H5::DataSet, double>(tra1DS, "atom_masses");
        std::vector<double> expectedAtomMasses = { 15.9989995956, 15.9989995956, 14.0066995621, 14.0066995621, 0.0 };
        CHECK_THAT(atomMasses, Catch::Matchers::Approx(expectedAtomMasses));

        auto nAtomsAttribute = readAttribute<H5::DataSet, std::size_t>(tra1DS, "number_of_atoms");
        std::vector<std::size_t> expectedAtomNumbers = {2, 3};
        CHECK_THAT(nAtomsAttribute, Catch::Matchers::Equals(expectedAtomNumbers));
    }
}

TEST_CASE("Test modularized force fields", "[CollisionModels][MDInteractionsModel]") {

    FileIO::MolecularStructureReader reader = FileIO::MolecularStructureReader();
    std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection = reader.readMolecularStructure("test_molecularstructure_reader.json");
    CollisionModel::Molecule ion({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, molecularStructureCollection.at("Ar+"));
    CollisionModel::Molecule background({3.0e-9, 0.0, 0.0}, {0.0, 0.0, 0.0}, molecularStructureCollection.at("Ar+"));

    std::vector<CollisionModel::Molecule*> molecules = {&ion, &background};
    std::vector<Core::Vector> forces = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    std::vector<Core::Vector> torque = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

    CollisionModel::MDForceField_LJ12_6 ff_lj_12_6 = CollisionModel::MDForceField_LJ12_6(0.208e-30);
    ff_lj_12_6.calculateForceField(molecules, forces, torque);

    CHECK(Approx(forces[0].x()).margin(1e-20) == 4.22513e-16);
    CHECK(Approx(forces[1].x()).margin(1e-20) == -4.22513e-16);
}


TEST_CASE("Test MD Interactions with rotation", "[CollisionModels][MDInteractionsModel]") {

    Core::globalRandomGeneratorPool = std::make_unique<Core::XoshiroTestRandomGeneratorPool>();

    FileIO::MolecularStructureReader reader = FileIO::MolecularStructureReader();
    std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection = reader.readMolecularStructure("test_molecularstructure_reader.json");
    Core::Particle ion;
    ion.setMolecularStructure(molecularStructureCollection.at("Ar+"));
    ion.setVelocity(Core::Vector(600.0, 50.0, 0.0));
    CollisionModel::MDForceField_LJ12_6 forceField(1.705E-30, "ALL", true);
    auto forceFieldPtr = std::make_unique<CollisionModel::MDForceField_LJ12_6>(forceField);
    CollisionModel::MDInteractionsModel mdSim = CollisionModel::MDInteractionsModel(2000000, 298,
                                                                                    28,
                                                                                    CollisionModel::MDInteractionsModel::DIAMETER_N2,
                                                                                    "N2",
                                                                                    1e-10, 
                                                                                    1E-17,
                                                                                    3, 1,
                                                                                    45e-10,
                                                                                    true,
                                                                                    std::move(forceFieldPtr),
                                                                                    molecularStructureCollection);

    double dt = 2e-11;
    std::string h5Filename = "MD_collisions_multiatom_trajectories_rotation_N2.h5";
    mdSim.setHDF5TrajectoryWriter(h5Filename, 45e-10, 0);
    mdSim.modifyVelocity(ion, dt);

    CHECK(Approx(ion.getVelocity().x()).margin(0.2) ==  495.7040847218);
    CHECK(Approx(ion.getVelocity().y()).margin(0.2) ==  43.9704656224);
    CHECK(Approx(ion.getVelocity().z()).margin(0.2) ==  -54.8414059866);

    unsigned int timestep = 0;
    double time = 0.0;
    for(int i = 0; i < 4; i++) {
        mdSim.updateModelTimestepParameters(timestep, time);
        mdSim.modifyVelocity(ion, 2e-11);
    }

    CHECK(Approx(ion.getVelocity().x()).margin(0.8) ==  32.3509170218);
    CHECK(Approx(ion.getVelocity().y()).margin(0.8) ==  -336.6975108956);
    CHECK(Approx(ion.getVelocity().z()).margin(0.2) ==  -326.1171405377);
}