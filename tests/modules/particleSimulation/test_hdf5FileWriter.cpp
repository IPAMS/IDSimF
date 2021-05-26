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
 test_hdf5FileWriter.cpp

 Testing of HDF5 trajectory file writer

 ****************************/

#include "PSim_trajectoryHDF5Writer.hpp"
#include "PSim_trajectoryExplorerJSONwriter.hpp"
#include "PSim_particleStartSplatTracker.hpp"
#include "Core_vector.hpp"
#include "BTree_particle.hpp"
#include "H5Cpp.h"
#include <array>
#include <vector>
#include <string>
#include <algorithm>
#include "catch.hpp"


// Define some helping data structures and functions .........
template <hsize_t NDIMS, typename DTYPE> struct DataField{
    hsize_t rank;
    std::array<hsize_t, NDIMS> dims;
    std::vector<DTYPE> data;

    DTYPE get(std::array<hsize_t,NDIMS> indices){
        hsize_t linIndex;
        if (NDIMS == 3) {
            linIndex = indices[0] * dims[2] * dims[1] + indices[1] * dims[2] + indices[2];
        }
        else if (NDIMS == 2){
            linIndex = indices[0] * dims[1] + indices[1];
        }
        else if(NDIMS == 1){
            linIndex = indices[0];
        }
        return data[linIndex];
    }
};

template <hsize_t NDIMS>DataField<NDIMS,double> readDataset(H5::DataSet& ds){

    //get the dataspace
    H5::DataSpace dataspace = ds.getSpace();

    //get the rank:
    int rank = dataspace.getSimpleExtentNdims();
    REQUIRE(NDIMS == rank);

    //get dimensions:
    hsize_t dims[NDIMS];
    int nDims = dataspace.getSimpleExtentDims(dims,NULL);
    REQUIRE(NDIMS == nDims);

    //prepare return object and prepare to read from HDF5 file:
    DataField<NDIMS,double> dField;
    dField.rank = nDims;
    hsize_t nElements = 1;
    hsize_t offset[nDims];
    //hsize_t count[nDims];

    for (hsize_t i=0; i<NDIMS; ++i){
        nElements *= dims[i];
        dField.dims[i] = dims[i];
        offset[i] = 0;
        //count[i] = dims[i];
    }

    //define selected hyperslab:
    dataspace.selectHyperslab(H5S_SELECT_SET, dims,offset);

    //define memory dataspace and hyperslab:
    H5::DataSpace memspace(NDIMS, dims);
    memspace.selectHyperslab(H5S_SELECT_SET,dims,offset);

    //read:
    double datBuf[nElements];
    const H5::PredType* nativeType;
    nativeType = &H5::PredType::NATIVE_DOUBLE;
    ds.read(datBuf,*nativeType,memspace,dataspace);

    for (hsize_t i=0; i<nElements; ++i){
        dField.data.emplace_back(datBuf[i]);
    }

    return dField;
}

std::vector<std::string> readStringAttribute(H5::Group& group, std::string attrName){
    H5::Attribute attr(group.openAttribute(attrName.c_str()));
    H5::DataSpace dataspace = attr.getSpace();

    //get dimensions:
    hsize_t dims[1];
    int nDims = dataspace.getSimpleExtentDims(dims, NULL);
    REQUIRE(nDims == 1);
    std::vector<std::string> result;
    char **datBuf = new char*[dims[0]];
    H5::StrType strdatatype(H5::PredType::C_S1, H5T_VARIABLE); // of length 256 characters
    attr.read(strdatatype, datBuf);
    for (hsize_t i = 0; i < dims[0]; ++i) {
        result.emplace_back(datBuf[i]);
    }
    return (result);
}

std::vector<int> readIntAttribute(H5::Group& group,std::string attrName){
    H5::Attribute attr(group.openAttribute(attrName.c_str()));
    H5::DataSpace dataspace = attr.getSpace();

    //get dimensions:
    hsize_t dims[1];
    int nDims = dataspace.getSimpleExtentDims(dims, NULL);
    REQUIRE(nDims == 1);
    std::vector<int> result;
    int datBuf[dims[0]];
    attr.read(H5::PredType::NATIVE_INT,datBuf);
    for (hsize_t i=0; i<dims[0]; ++i){
        result.emplace_back(datBuf[i]);
    }
    return result;
}

std::vector<double> readDoubleAttribute(H5::Group& group,std::string attrName){
    H5::Attribute attr(group.openAttribute(attrName.c_str()));
    H5::DataSpace dataspace = attr.getSpace();

    //get dimensions:
    hsize_t dims[1];
    int nDims = dataspace.getSimpleExtentDims(dims, NULL);
    REQUIRE(nDims == 1);
    std::vector<double> result;
    double datBuf[dims[0]];
    attr.read(H5::PredType::NATIVE_DOUBLE,datBuf);
    for (hsize_t i=0; i<dims[0]; ++i){
        result.emplace_back(datBuf[i]);
    }
    return result;
}


TEST_CASE( "Test HDF5 trajectory file writer", "[ParticleSimulation][file writers]") {
    std::string filenameBare("test_particles_bare.h5");
    std::string filenameAux("test_particles_aux.h5");
    int nParticles = 5;
    int nFrames = 8;

    SECTION( "Hdf5 trajectory writer can write trajectory files with and without particle attributes"){

        // prepare and write data to a hdf5 file:
        //prepare a bare writer:
        ParticleSimulation::TrajectoryHDF5Writer writerBare(filenameBare, false);

        //prepare a hdf5 writer including parameter attributes
        ParticleSimulation::partAttribTransformFctType pAttribTransformFct =
                [](BTree::Particle *particle) -> std::vector<double>{
                    std::vector<double> result = {
                            particle->getVelocity().x(),
                            particle->getVelocity().y(),
                            particle->getVelocity().z(),
                            particle->getAcceleration().x(),
                            particle->getAcceleration().y(),
                            particle->getAcceleration().z()
                    };
                    return result;
                };
        std::vector<std::string> pAttribNames = {"velocity x", "velocity y", "velocity z",
                                                 "acceleration x", "acceleration y", "acceleration z"};

        //prepare a hdf5 writer including parameter attributes
        ParticleSimulation::partAttribTransformFctTypeInteger pAttribTransformFctInt =
                [](BTree::Particle *particle) -> std::vector<int>{
                    std::vector<int> result = {
                            particle->getIntegerAttribute("global index")
                    };
                    return result;
                };
        std::vector<std::string> pAttribNamesInt = {"global index"};

        ParticleSimulation::TrajectoryHDF5Writer writerAux(filenameAux);
        writerAux.setParticleAttributes(pAttribNames, pAttribTransformFct);
        writerAux.setParticleAttributes(pAttribNamesInt, pAttribTransformFctInt);

        // prepare data structures:
        std::vector<BTree::uniquePartPtr>particles;
        std::vector<BTree::Particle *> particlePtrs;

        // write an empty frame:
        writerBare.writeTimestep(particlePtrs, 0.0);
        writerAux.writeTimestep(particlePtrs, 0.0);

        //prepare particles to test:
        ParticleSimulation::ParticleStartSplatTracker tracker;
        for (int i=0; i<nParticles; ++i){
            double timeOfBirth = i*0.01;
            BTree::uniquePartPtr particle = std::make_unique<BTree::Particle>();
            particle->setLocation({0.0, 1.0, 0.0});
            particlePtrs.emplace_back(particle.get());
            particles.emplace_back(std::move(particle));
            tracker.particleStart(particlePtrs.at(i), timeOfBirth);
        }

        for (int k = 1; k < nFrames; k++) {
            for (int i = 0; i < nParticles; ++i) {
                particles[i]->setLocation(Core::Vector(i * 0.1, 1.0, k * 10.0));
                particles[i]->setAcceleration(Core::Vector(i * 0.1, k * 1.0, 0));
                particles[i]->setVelocity(Core::Vector(i * 0.01, i * 0.1, k * 10.0));
            }
            writerBare.writeTimestep(particlePtrs, k * 1.0);
            writerAux.writeTimestep(particlePtrs, k * 1.0);
        }
        tracker.particleSplat(particlePtrs[0], 1.0);
        tracker.particleSplat(particlePtrs[1], nFrames*1.0);
        tracker.particleSplat(particlePtrs[nParticles-1], 0.99);

        writerBare.writeStartSplatData(tracker);
        writerBare.finalizeTrajectory();

        std::vector<double> additionalVectorDataset = {10.0,10.1,10.2,10.3,10.4};
        writerAux.writeNumericListDataset("additional_vector_data_set", additionalVectorDataset);

        std::vector<std::string> particleNames = {"Substance 1","Substance 2"};
        writerAux.writeTrajectoryAttribute("particle names",particleNames);

        std::vector<double> particleMasses = {125.5,128.8};
        writerAux.writeTrajectoryAttribute("particle masses",particleMasses);
        writerAux.finalizeTrajectory();
    }

    SECTION("Written bare hdf5 trajectory contains correct data"){
        //read that file and check contents:
        H5::H5File bareFile(filenameBare.c_str(),H5F_ACC_RDONLY);

        //check file version:
        H5::Group group (bareFile.openGroup("particle_trajectory"));
        auto fileVersionId = readIntAttribute(group,"file version");
        REQUIRE_NOTHROW(fileVersionId[0] == 2);

        //check timesteps:
        H5::DataSet dsTimes= bareFile.openDataSet("particle_trajectory/times");
        auto dFieldTimes = readDataset<1>(dsTimes);
        REQUIRE(dFieldTimes.rank == 1);
        std::vector<double> times = dFieldTimes.data;
        int nTimesteps = times.size();

        REQUIRE(nTimesteps == nFrames);

        std::array<hsize_t,1> index= {0};
        for (int ts=1; ts<nTimesteps; ++ts) {
            index[0] = ts;
            REQUIRE(Approx(dFieldTimes.get(index)) == ts * 1.0);
            std::string tsPath = "/particle_trajectory/timesteps/" + std::to_string(ts) +"/positions";

            H5::DataSet dsPositions = bareFile.openDataSet(tsPath.c_str());
            auto dField = readDataset<2>(dsPositions);
            REQUIRE(dField.rank == 2);

            hsize_t nParticles = dField.dims[0];
            hsize_t nSpatialDims = dField.dims[1];

            REQUIRE(nParticles == 5);
            REQUIRE(nSpatialDims == 3);

            //check position data:
            std::array<hsize_t,2> indices = {0,0};
            for (hsize_t pi = 0; pi < nParticles; ++pi) {
                indices[0] = pi;
                indices[1] = 0;
                REQUIRE(Approx(dField.get(indices)) == 0.1 * pi);
                indices[1] = 1;
                REQUIRE(Approx(dField.get(indices)) == 1.0);
                indices[1] = 2;
                REQUIRE(Approx(dField.get(indices)) == ts*10);
            }
        }

        //check splattimes:
        H5::DataSet dsStartTimes= bareFile.openDataSet("particle_trajectory/start_splat/particle start times");
        H5::DataSet dsSplatTimes= bareFile.openDataSet("particle_trajectory/start_splat/particle splat times");

        auto dFieldSplattimes= readDataset<2>(dsSplatTimes);
        REQUIRE(dFieldSplattimes.rank == 2);
        std::array<hsize_t,2> indices = {0,0};
        REQUIRE(Approx(dFieldSplattimes.get(indices)) == 1.0);
        indices[0] = 1;
        REQUIRE(Approx(dFieldSplattimes.get(indices)) == nFrames*1.0);
        indices[0] = 2;
        REQUIRE(Approx(dFieldSplattimes.get(indices)) == 0.0);
        indices[0] = nParticles-1;
        REQUIRE(Approx(dFieldSplattimes.get(indices)) == 0.99);
    }

    SECTION("Written hdf5 trajectory contains correct particle attributes"){
        H5::H5File auxFile(filenameAux.c_str(),H5F_ACC_RDONLY);
        H5::DataSet dsTimes= auxFile.openDataSet("particle_trajectory/times");
        auto dFieldTimes = readDataset<1>(dsTimes);
        REQUIRE(dFieldTimes.rank == 1);
        std::vector<double> times = dFieldTimes.data;
        int nTimesteps = times.size();

        REQUIRE(nTimesteps == nFrames);

        //read that file and check contents:
        for (int ts=1; ts<nTimesteps; ++ts) {

            //check float particle attributes:
            std::string pAttribFloatPath = "/particle_trajectory/timesteps/" + std::to_string(ts) +"/particle_attributes_float";

            H5::DataSet dsAttribFloat= auxFile.openDataSet(pAttribFloatPath.c_str());
            auto dField = readDataset<2>(dsAttribFloat);
            REQUIRE(dField.rank == 2);

            hsize_t nParticles = dField.dims[0];
            hsize_t nAttributeDims = dField.dims[1];

            REQUIRE(nParticles == 5);
            REQUIRE(nAttributeDims == 6);
            std::array<hsize_t,2> indices = {0,0};

            for (hsize_t pi = 0; pi < nParticles; ++pi) {
                indices[0] = pi;
                indices[1] = 0;
                REQUIRE(Approx(dField.get(indices)) == 0.01 * pi);
                indices[1] = 4;
                REQUIRE(Approx(dField.get(indices)) == 1.0 * ts);
                indices[1] = 5;
                REQUIRE(Approx(dField.get(indices)) == 0.0);
            }


            //check integer particle attributes:
            std::string pAttribIntegerPath = "/particle_trajectory/timesteps/" + std::to_string(ts) +"/particle_attributes_integer";

            H5::DataSet dsAttribInteger= auxFile.openDataSet(pAttribIntegerPath.c_str());
            auto dFieldInteger = readDataset<2>(dsAttribInteger);
            REQUIRE(dFieldInteger.rank == 2);

            nParticles = dFieldInteger.dims[0];
            nAttributeDims = dFieldInteger.dims[1];

            REQUIRE(nParticles == 5);
            REQUIRE(nAttributeDims == 1);
            indices = {0,0};

            for (hsize_t pi = 0; pi < nParticles; ++pi) {
                indices[0] = pi;
                REQUIRE(dFieldInteger.get(indices) == pi);
            }
        }

        H5::DataSet dsAdditional= auxFile.openDataSet("particle_trajectory/additional_vector_data_set");
        auto dFieldAdditional = readDataset<2>(dsAdditional);
        REQUIRE(dFieldAdditional.rank == 2);
        hsize_t nAdditional = dFieldAdditional.dims[0];
        REQUIRE(nAdditional == 5);

        //check additional data:
        std::array<hsize_t,2> indicesAdditional = {0,0};
        for (hsize_t i=0; i<nAdditional; ++i) {
            indicesAdditional[0] = i;
            REQUIRE(Approx(dFieldAdditional.get(indicesAdditional)) == 10 + i*0.1);
        }

        //check additional attributes:
        H5::Group group (auxFile.openGroup("particle_trajectory"));

        auto attribParticleNames = readStringAttribute(group,"particle names");
        REQUIRE(attribParticleNames.size() == 2);
        REQUIRE(attribParticleNames[0] == "Substance 1");
        REQUIRE(attribParticleNames[1] == "Substance 2");

        auto attribParticleMasses = readDoubleAttribute(group,"particle masses");
        REQUIRE(attribParticleMasses.size() == 2);
        REQUIRE(Approx(attribParticleMasses[0]) == 125.5);
        REQUIRE(Approx(attribParticleMasses[1]) == 128.8);
    }
}