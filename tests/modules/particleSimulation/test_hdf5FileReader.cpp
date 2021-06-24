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
 test_hdf5FileReader.cpp

 Testing of HDF5 Reader / Wrapper class

 ****************************/


#include "PSim_HDF5Reader.hpp"
#include <vector>
#include <array>
#include "catch.hpp"

TEST_CASE("Test HDF5 file reader", "[ParticleSimulation][file reader]") {

    ParticleSimulation::HDF5Reader h5reader_trajectory("test_particle_trajectory.h5");
    ParticleSimulation::HDF5Reader h5reader_scalarField("test_linear_scalar_field_01.h5");

    SECTION("Testing of attribute reading") {
        SECTION("Hdf5 reader can read string attribute vector") {
            std::vector<std::string> pNames = h5reader_trajectory.readAttributeVector<std::string>(
                    "particle_trajectory", "particle names");
            REQUIRE(pNames.at(0)=="Substance 1");
            REQUIRE(pNames.at(1)=="Substance 2");
            REQUIRE(pNames.size()==2);
        }

        SECTION("Hdf5 reader can read integer attribute vector") {
            std::vector<int> timeSteps = h5reader_trajectory.readAttributeVector<int>(
                    "particle_trajectory", "number of timesteps");
            REQUIRE(timeSteps.size()==1);
            REQUIRE(timeSteps[0]==8);
        }

        SECTION("Hdf5 reader can read double attribute vector") {
            std::vector<double> pMasses = h5reader_trajectory.readAttributeVector<double>(
                    "particle_trajectory", "particle masses");
            REQUIRE(pMasses.size()==2);
            REQUIRE(Approx(pMasses[0])==125.5);
            REQUIRE(Approx(pMasses[1])==128.8);
        }

        SECTION("Reading of attribute vector with wrong datatype throws exception") {
            REQUIRE_THROWS(h5reader_trajectory.readAttributeVector<int>("particle_trajectory", "particle names"));
            REQUIRE_THROWS(h5reader_trajectory.readAttributeVector<std::string>(
                    "particle_trajectory", "particle masses"));
        }
    }

    SECTION("Testing of metadata reading") {
        SECTION("Reading of object number in group works") {
            hsize_t numOfObjects = h5reader_trajectory.numberOfObjectsInGroup("particle_trajectory/timesteps");
            REQUIRE(numOfObjects==8);
        }

        SECTION("Reading of dataset dimensionality works") {
            int nDimsTraj = h5reader_trajectory.datasetNDims("particle_trajectory/timesteps/0/aux_parameters");
            REQUIRE(nDimsTraj==2);

            int nDimsField = h5reader_scalarField.datasetNDims("fields/test_field");
            REQUIRE(nDimsField==3);
        }

        SECTION("Reading of object and dataset names in group works") {
            std::vector<std::string> objectNames = h5reader_trajectory.namesOfObjectsInGroup(
                    "particle_trajectory/timesteps");
            std::vector<std::string> correctObjectNames = {"0", "1", "2", "3", "4", "5", "6", "7"};
            REQUIRE(objectNames==correctObjectNames);

            std::vector<std::string> datasetNames = h5reader_trajectory.namesOfDatasetsInGroup(
                    "particle_trajectory/timesteps");
            REQUIRE(datasetNames.empty());

            std::vector<std::string> datasetNamesInFrame = h5reader_trajectory.namesOfDatasetsInGroup(
                    "particle_trajectory/timesteps/"+objectNames[0]);
            std::vector<std::string> correctDatasetNames = {"aux_parameters", "positions"};
            REQUIRE(datasetNamesInFrame==correctDatasetNames);
        }
    }

    SECTION("Data sets with double type are readable"){

        SECTION("Read times vector from particle trajectory"){
            ParticleSimulation::HDF5Reader::DataField dFieldTimes = h5reader_trajectory.readDataset<1>(
                    "particle_trajectory/times");
            REQUIRE(dFieldTimes.rank == 1);
            REQUIRE(dFieldTimes.data.size() == 8);
            std::array<hsize_t,1> indices = {6};
            REQUIRE(dFieldTimes.get(indices) == Approx(6.0));
        }

        SECTION("Read particle simulation frame from particle trajectory"){
            ParticleSimulation::HDF5Reader::DataField auxParams =
                    h5reader_trajectory.readDataset<2>("particle_trajectory/timesteps/3/aux_parameters");
            REQUIRE(auxParams.rank == 2);
            REQUIRE(auxParams.get({2,1}) == Approx(0.2));

            ParticleSimulation::HDF5Reader::DataField positions =
                    h5reader_trajectory.readDataset<2>("particle_trajectory/timesteps/3/positions");
            REQUIRE(positions.rank == 2);
            REQUIRE(positions.get({4,0}) == Approx(0.4));
        }

        SECTION("Read scalar field data set"){
            ParticleSimulation::HDF5Reader::DataField scalarField =
                    h5reader_scalarField.readDataset<3>("fields/test_field");
            REQUIRE(scalarField.rank == 3);
            REQUIRE(scalarField.get({1,2,3}) == Approx(19.0));
            REQUIRE(scalarField.get({0,2,1}) == Approx(12.0));
        }

        SECTION("Read small vector field data set"){
            ParticleSimulation::HDF5Reader h5reader_vectorField("test_linear_vector_field_01.h5");

            ParticleSimulation::HDF5Reader::DataField vecField_1 =
                    h5reader_vectorField.readDataset<4>("fields/test_vectorfield_1");

            ParticleSimulation::HDF5Reader::DataField vecField_2 =
                    h5reader_vectorField.readDataset<4>("fields/test_vectorfield_2");

            REQUIRE(vecField_1.rank == 4);
            REQUIRE(vecField_1.dims == std::array<hsize_t, 4>({12, 3, 2, 3}));

            REQUIRE(vecField_2.rank == 4);
            REQUIRE(vecField_2.dims == std::array<hsize_t, 4>({12, 3, 2, 3}));

            REQUIRE(vecField_1.get({2,1,0,0}) == Approx(2));
            REQUIRE(vecField_1.get({2,1,0,1}) == Approx(5));
            REQUIRE(vecField_1.get({2,1,0,2}) == Approx(1));

            REQUIRE(vecField_2.get({0,2,0,0}) == Approx(10));
            REQUIRE(vecField_2.get({0,2,0,1}) == Approx(15));
            REQUIRE(vecField_2.get({0,2,1,2}) == Approx(11));
            REQUIRE(vecField_2.get({0,0,0,0}) == Approx(-10));
        }

        SECTION("Read large vector field file") {
            ParticleSimulation::HDF5Reader h5reader_flow("quad_dev_flow_3d.h5");
            auto ds = h5reader_flow.readDataset<4>("/fields/velocity");
            REQUIRE(ds.dims == std::array<hsize_t, 4>({300, 40, 40, 3}));
            REQUIRE(ds.data.size() == 300* 40* 40* 3);
        }

        //TODO: test dataset reading with wrong rank
        //TODO: test dataset reading with non existing dataset
    }
}