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
 ****************************/

#include "PSim_trajectoryHDF5Writer.hpp"
#include "PSim_particleStartSplatTracker.hpp"
#include <array>
#include <cmath>
#include <string>

/**
 * Constructs a new HDF5 trajectory filewriter
 * @param hdf5Filename A filename / path of a HDF5 file to write
 * @param compression If true: the HDF5 file is written with data compression
 */
ParticleSimulation::TrajectoryHDF5Writer::TrajectoryHDF5Writer(const std::string& hdf5Filename, bool compression):
        compression_(compression),
        memspaceTimestep_(1, slabDimsTimestep_)
{
    h5f_ = std::make_unique<H5::H5File>(hdf5Filename.c_str(), H5F_ACC_TRUNC);
    baseGroup_ = std::make_unique<H5::Group>(h5f_->createGroup("particle_trajectory"));
    h5f_->createGroup("/particle_trajectory/timesteps");

    //write version number
    writeTrajectoryAttribute("file version",FILE_TYPE_VERSION);

    //prepare timestep dataset structures:
    hsize_t dimsTimesteps[1] = {0};                       // dataset dimensions at creation
    hsize_t maxdimsTimesteps[1] = {H5S_UNLIMITED};        // maximum dataset dimensions
    hsize_t chunkDimsTimesteps[1] = {5};
    H5::DataSpace dataspaceTimesteps(1,dimsTimesteps,maxdimsTimesteps);

    // Modify dataset creation properties to enable chunking and optional compression:
    H5::DSetCreatPropList propTimesteps;
    propTimesteps.setChunk(1, chunkDimsTimesteps);

    //create actual dataset for times:
    dsetTimesteps_ = std::make_unique<H5::DataSet>(h5f_->createDataSet("/particle_trajectory/times",
                                                                H5::PredType::IEEE_F32BE, dataspaceTimesteps,
                                                                propTimesteps));
}


/**
 * FIXME
 *
 * @param attributeNames
 * @param attributesTransformFct
 */
void ParticleSimulation::TrajectoryHDF5Writer::setParticleAttributes(const std::vector<std::string>& attributeNames,
                                                                     partAttribTransformFctType attributesTransformFct) {
    hasParticleAttributes_ = true;
    particleAttributeTransformFct_ = std::move(attributesTransformFct);
    nPAttributes_ = attributeNames.size();
    writeTrajectoryAttribute("attributes names",attributeNames);
}

/**
 *
 * @param attributeNames
 * @param attributesTransformFct
 */
void ParticleSimulation::TrajectoryHDF5Writer::setParticleAttributes(const std::vector<std::string>& attributeNames,
                                                                     partAttribTransformFctTypeInteger attributesTransformFct) {
    hasParticleAttributesInteger_ = true;
    particleAttributeTransformFctInteger_ = std::move(attributesTransformFct);
    nPAttributesInteger_ = attributeNames.size();
    writeTrajectoryAttribute("integer attributes names",attributeNames);
}

/**
 * Writes a single time step to the HDF5 trajectory file
 *
 * Note that an empty timestep, without any particles, is **silently ignored**.
 *
 * @param particles The particle ensemble in the simulation to write to the trajectory file
 * @param time The current simulated time
 */
void ParticleSimulation::TrajectoryHDF5Writer::writeTimestep(std::vector<BTree::Particle*> &particles,
                                                             double time){

    // Write time of this time step to the times vector:
    sizeTimesteps_[0] += 1;
    dsetTimesteps_->extend(sizeTimesteps_);
    H5::DataSpace dSpaceTimestep = dsetTimesteps_->getSpace();
    dSpaceTimestep.selectHyperslab(H5S_SELECT_SET, slabDimsTimestep_, offsetScalarLike_);
    double ts[] = {time};
    dsetTimesteps_->write(ts, H5::PredType::NATIVE_DOUBLE, memspaceTimestep_, dSpaceTimestep);

    // Create group for this timestep and dataset for the location data:
    std::string timeStepGroupPath = "/particle_trajectory/timesteps/" + std::to_string(sizeTimesteps_[0]-1);
    timeStepGroup_ = std::make_unique<H5::Group>(h5f_->createGroup(timeStepGroupPath.c_str()));


    hsize_t nParticles = particles.size();
    // Write particle number as attribute to the time step:
    writeAttribute_(timeStepGroup_, "number of particles", nParticles);

    if (nParticles > 0) {
        //write particle location data:

        //define chunk size in parameter direction:
        hsize_t nParticlesChunk = 256;
        if (nParticlesChunk>nParticles) {
            nParticlesChunk = nParticles;
        }

        //prepare location dataset structures:
        hsize_t dimsLocation[2] = {nParticles, 3};            // dataset dimensions at creation
        hsize_t maxdimsLocation[2] = {nParticles, 3};         // maximum dataset dimensions
        hsize_t chunkDimsLocation[2] = {nParticlesChunk, 3};
        H5::DataSpace dataspaceLocation(2, dimsLocation, maxdimsLocation);

        // Modify dataset creation properties to enable chunking and optional compression:
        H5::DSetCreatPropList propLocation;
        propLocation.setChunk(2, chunkDimsLocation);

        if (compression_) {
            propLocation.setDeflate(6);
        }

        // Create dataset for the location data:
        std::unique_ptr<H5::DataSet> dsetPPositions = std::make_unique<H5::DataSet>(
                timeStepGroup_->createDataSet("positions",
                        H5::PredType::IEEE_F32BE, dataspaceLocation,
                        propLocation));

        //prepare dataset:
        H5::DataSpace dSpaceLocation = dsetPPositions->getSpace();
        hsize_t slabDimsLocation[2] = {nParticles, 3};
        hsize_t offsetLocation[2] = {0, 0};
        dSpaceLocation.selectHyperslab(H5S_SELECT_SET, slabDimsLocation, offsetLocation);

        //create location data buffer:
        std::vector<double> bufLocation(nParticles*3);
        for (hsize_t i = 0; i<nParticles; ++i) {
            Core::Vector loc = particles[i]->getLocation();
            bufLocation[i*3] = loc.x();
            bufLocation[i*3+1] = loc.y();
            bufLocation[i*3+2] = loc.z();
        }

        //write to the dataset:
        H5::DataSpace memspaceLocation(2, slabDimsLocation);
        dsetPPositions->write(bufLocation.data(), H5::PredType::NATIVE_DOUBLE, memspaceLocation, dSpaceLocation);

        if (hasParticleAttributes_) {
            writeTimestepParticleAttributes_(particles);
        }
        if (hasParticleAttributesInteger_) {
            writeTimestepParticleAttributesInteger_(particles);
        }
    }
    offsetScalarLike_[0] +=1;
}

/**
 * Write a numeric vector to a dataset in the trajectory file
 * @param dsName Name of the dataset in the HDF5 file to write the vector into
 * @param values Numeric vector with values to write to the dataset
 * @param group FIXME
 */

template <typename DT>
void ParticleSimulation::TrajectoryHDF5Writer::writeNumericListDataset(std::string dsName, const std::vector<DT> &values, H5::Group* group){
    std::vector<std::array<DT, 1>> valuesPacked;
    for (auto const &val: values){
        std::array<DT,1> ar = {val};
        valuesPacked.emplace_back(ar);
    }
    writeArrayDataSet<DT, 1>(dsName, valuesPacked, group);
};

void ParticleSimulation::TrajectoryHDF5Writer::write3DVectorListDataset(std::string dsName, const std::vector<Core::Vector> &values, H5::Group* group){
    std::vector<std::array<double, 3>> valuesPacked;

    for (auto const &val: values){
        std::array<double,3> ar = {val.x(), val.y(), val.z()};
        valuesPacked.emplace_back(ar);
    }
    writeArrayDataSet<double, 3>(dsName, valuesPacked, group);
};


/**
 * FIXME
 * @tparam DT
 * @tparam NCOLUMNS
 * @param dsName
 * @param values
 * @param group
 */
template <typename DT, int NCOLUMNS>
void ParticleSimulation::TrajectoryHDF5Writer::writeArrayDataSet(std::string dsName, const std::vector<std::array<DT, NCOLUMNS>> &values, H5::Group* group){

    //prepare dataset structures:
    size_t nValues = values.size();

    hsize_t dims[2] = {nValues, NCOLUMNS};                       // dataset dimensions at creation
    hsize_t offset[2] = {0, 0};
    //hsize_t maxdims[2] = {nValues, NCOLUMNS};                   // maximum dataset dimensions gives
    hsize_t chunkDims[2] = {nValues, NCOLUMNS};
    H5::DataSpace dataspace(2, dims);

    // Modify dataset creation properties to enable chunking and optional compression:
    H5::DSetCreatPropList props;
    props.setChunk(2, chunkDims);

    //create actual datasets:
    if (group ==nullptr){
        group = baseGroup_.get();
    }

    const H5::PredType* dsetPredType;
    const H5::PredType* writePredType;
    if constexpr (std::is_same<DT, double>::value) {
        dsetPredType = &H5::PredType::IEEE_F32BE;
        writePredType = &H5::PredType::NATIVE_DOUBLE;
    }
    else if constexpr (std::is_same<DT, int>::value) {
        dsetPredType = &H5::PredType::NATIVE_INT;
        writePredType = &H5::PredType::NATIVE_INT;
    }

    H5::DataSet dset = group->createDataSet(dsName.c_str(), *dsetPredType, dataspace, props);
    H5::DataSpace dspace = dset.getSpace();
    dspace.selectHyperslab(H5S_SELECT_SET, dims, offset);

    // Define memory space.
    std::vector<DT> datBuf(nValues*NCOLUMNS);
    for (std::size_t i=0; i<nValues; ++i){
        for(int j=0; j<NCOLUMNS; ++j){
            datBuf[i*NCOLUMNS+j] = values[i][j];
        }
    }
    H5::DataSpace memspace(2, dims);
    dset.write(datBuf.data(), *writePredType, memspace, dataspace);
}


/**
 * Writes an integer attribute to the trajectory file
 * @param attrName Name of the attribute to write
 * @param value The value to write into the attribute in the trajectory file
 */
void ParticleSimulation::TrajectoryHDF5Writer::writeTrajectoryAttribute(std::string attrName, int value){

    writeAttribute_(baseGroup_, attrName, value);
}

/**
 * Writes an numeric vector attribute to the trajectory file
 * @param attrName Name of the attribute to write
 * @param value The values to write into the attribute in the trajectory file
 */
void ParticleSimulation::TrajectoryHDF5Writer::writeTrajectoryAttribute(std::string attrName,
                                                                        const std::vector<double> &values){
    hsize_t nVals = values.size();
    hsize_t dims[1] = { nVals };
    H5::DataSpace attr_dataspace = H5::DataSpace (1, dims);
    H5::Attribute doubleAttribute = baseGroup_->createAttribute(attrName.c_str(), H5::PredType::IEEE_F32BE, attr_dataspace);

    double data[nVals];
    for (size_t i = 0; i < nVals; ++i)
    {
        data[i] = values[i];
    }
    doubleAttribute.write(H5::PredType::NATIVE_DOUBLE, data);
}

/**
 * Writes an string vector attribute to the trajectory file
 * @param attrName Name of the attribute to write
 * @param value The values to write into the attribute in the trajectory file
 */
void ParticleSimulation::TrajectoryHDF5Writer::writeTrajectoryAttribute(std::string attrName,
                                                                        const std::vector<std::string> &values){
    hsize_t nVals = values.size();
    hsize_t dims[1] = { nVals };
    H5::DataSpace attr_dataspace = H5::DataSpace (1, dims);

    // Create new string datatype for attribute
    H5::StrType strdatatype(H5::PredType::C_S1, 256); // of length 256 characters
    strdatatype.setSize(H5T_VARIABLE); //strings can be of variable size...

    // Create attribute to write to
    H5::Attribute stringAttribute = baseGroup_->createAttribute(attrName.c_str(), strdatatype, attr_dataspace);

    //Due to a abi bug between the std lib and hdf, we need to provide raw strings to the hdf5 methods
    //Save array of pointers to the raw c strings and use that array as data buffer
    const char* data[nVals];
    for (size_t i = 0; i < nVals; ++i)
    {
        data[i] = values[i].c_str();
    }
    stringAttribute.write(strdatatype, data);
}

/**
 * Writes the splat (termination) times of the simulated particle ensemble to the HDF5 trajectory file
 * @param particles The simulated particle ensemble
 */
void ParticleSimulation::TrajectoryHDF5Writer::writeSplatTimes(std::vector<BTree::Particle *> &particles){
    //prepare timestep dataset structures:
    hsize_t nParticles = particles.size();

    hsize_t dimsSplattimes[1] = {nParticles};                        // dataset dimensions at creation
    hsize_t maxdimsSplattimes[1] = {nParticles};                     // maximum dataset dimensions

    //define chunk size in parameter direction:
    hsize_t nParticlesChunk = 200;
    if (nParticlesChunk > nParticles){
        nParticlesChunk = nParticles;
    }

    hsize_t chunkDimsSplattimes[1] = {nParticlesChunk};              // chunk for writing
    H5::DataSpace dataspaceSplattimes(1,dimsSplattimes,maxdimsSplattimes);

    // Modify dataset creation properties to enable chunking and optional compression:
    H5::DSetCreatPropList propSplattimes;
    propSplattimes.setChunk(1, chunkDimsSplattimes);

    //create actual datasets:
    std::unique_ptr<H5::DataSet> dsetSplattimes = std::make_unique<H5::DataSet>(h5f_->createDataSet("/particle_trajectory/splattimes",
                                                                       H5::PredType::IEEE_F32BE, dataspaceSplattimes,
                                                                       propSplattimes));

    // Define memory space and write.
    hsize_t slabDimsSplattimes[1] = {nParticles};
    H5::DataSpace memspaceSplattimes(1, slabDimsSplattimes);

    std::vector<double> datBuf(nParticles); //use a vector to prevent stack overflows
    for(hsize_t i=0; i<nParticles; ++i){
        datBuf[i] = particles[i]->getSplatTime();
    }
    dsetSplattimes->write(datBuf.data(),H5::PredType::NATIVE_DOUBLE,memspaceSplattimes,dataspaceSplattimes);
}

void ParticleSimulation::TrajectoryHDF5Writer::writeStartSplatData(ParticleStartSplatTracker tracker) {

    tracker.sortStartSplatData();

    H5::Group startSplatGroup = baseGroup_->createGroup("start_splat");
    writeNumericListDataset("particle splat state", tracker.getSplatState(), &startSplatGroup);
    writeNumericListDataset("particle start times", tracker.getStartTimes(), &startSplatGroup);
    writeNumericListDataset("particle splat times", tracker.getSplatTimes(), &startSplatGroup);
    write3DVectorListDataset("particle start locations", tracker.getStartLocations(), &startSplatGroup);
    write3DVectorListDataset("particle splat locations", tracker.getSplatLocations(), &startSplatGroup);
}

/**
 * Finalizes the trajectory, usually after the simulation has finished
 */
void ParticleSimulation::TrajectoryHDF5Writer::finalizeTrajectory(){
    writeTrajectoryAttribute("number of timesteps",offsetScalarLike_[0]);
}

/**
 * Writes the additional particle attributes of a time step
 * @param particles The simulated particle ensemble to write the auxiliary data for
 */
void ParticleSimulation::TrajectoryHDF5Writer::writeTimestepParticleAttributes_(std::vector<BTree::Particle*> &particles){

    hsize_t nParticles = particles.size();

    //define chunk size in parameter direction:
    hsize_t nParticlesChunk = 200;
    if (nParticlesChunk > nParticles){
        nParticlesChunk = nParticles;
    }

    //prepare location dataset structures:
    hsize_t dimsAux[2] = {nParticles, nPAttributes_};            // dataset dimensions at creation
    hsize_t maxdimsAux[2] = {nParticles, nPAttributes_};         // maximum dataset dimensions
    hsize_t chunkDimsAux[2] = {nParticlesChunk, nPAttributes_};
    H5::DataSpace dataspaceAux(2,dimsAux,maxdimsAux);

    H5::DSetCreatPropList propAux;
    propAux.setChunk(2, chunkDimsAux);

    //create and prepare dataset for aux parameters:
    std::unique_ptr<H5::DataSet> dset = std::make_unique<H5::DataSet>(
            timeStepGroup_->createDataSet("particle_attributes_float",
                    H5::PredType::IEEE_F32BE, dataspaceAux,
                    propAux));

    H5::DataSpace dSpaceAttrib = dset ->getSpace();
    hsize_t slabDimsLocation[2]  = {nParticles, nPAttributes_};
    hsize_t offsetLocation[2]  = {0,0};
    dSpaceAttrib.selectHyperslab(H5S_SELECT_SET, slabDimsLocation, offsetLocation);

    //create and fill aux data buffer:
    std::vector<double> bufAux(nParticles*nPAttributes_);
    for (hsize_t i=0; i<nParticles; ++i){
        std::vector<double> auxDat = particleAttributeTransformFct_(particles[i]);
        for(hsize_t j=0; j<nPAttributes_; ++j){
            bufAux[i*nPAttributes_+j] = auxDat[j];
        }
    }

    //write to the dataset:
    H5::DataSpace memspaceAux(2,slabDimsLocation);
    dset->write(bufAux.data(),H5::PredType::NATIVE_DOUBLE,memspaceAux,dSpaceAttrib);
}

/**
 * FIXME
 * @param particles The simulated particle ensemble to write the auxiliary data for
 */
void ParticleSimulation::TrajectoryHDF5Writer::writeTimestepParticleAttributesInteger_(std::vector<BTree::Particle*> &particles){

    hsize_t nParticles = particles.size();

    //define chunk size in parameter direction:
    hsize_t nParticlesChunk = 200;
    if (nParticlesChunk > nParticles){
        nParticlesChunk = nParticles;
    }

    //prepare location dataset structures:
    hsize_t dimsAux[2] = {nParticles, nPAttributesInteger_};            // dataset dimensions at creation
    hsize_t maxdimsAux[2] = {nParticles, nPAttributesInteger_};         // maximum dataset dimensions
    hsize_t chunkDimsAux[2] = {nParticlesChunk, nPAttributesInteger_};
    H5::DataSpace dataspaceAux(2,dimsAux,maxdimsAux);

    H5::DSetCreatPropList propAux;
    propAux.setChunk(2, chunkDimsAux);

    //create and prepare dataset for aux parameters:
    std::unique_ptr<H5::DataSet> dset = std::make_unique<H5::DataSet>(
            timeStepGroup_->createDataSet("particle_attributes_integer",
                    H5::PredType::NATIVE_INT, dataspaceAux,
                    propAux));

    H5::DataSpace dSpaceAttrib = dset ->getSpace();
    hsize_t slabDimsLocation[2]  = {nParticles, nPAttributesInteger_};
    hsize_t offsetLocation[2]  = {0,0};
    dSpaceAttrib.selectHyperslab(H5S_SELECT_SET, slabDimsLocation, offsetLocation);

    //create and fill aux data buffer:
    std::vector<int> bufAux(nParticles*nPAttributesInteger_);
    for (hsize_t i=0; i<nParticles; ++i){
        std::vector<int> auxDat = particleAttributeTransformFctInteger_(particles[i]);
        for(hsize_t j=0; j<nPAttributesInteger_; ++j){
            bufAux[i*nPAttributesInteger_+j] = auxDat[j];
        }
    }

    //write to the dataset:
    H5::DataSpace memspaceAux(2, slabDimsLocation);
    dset->write(bufAux.data(),H5::PredType::NATIVE_INT, memspaceAux, dSpaceAttrib);
}


/**
 * Writes an integer attribute to a HDF5 group
 * @param group The group to write to
 * @param attrName Name of the attribute to write
 * @param value The value to write into the attribute in the trajectory file
 */
void ParticleSimulation::TrajectoryHDF5Writer::writeAttribute_(std::unique_ptr<H5::Group>& group, const std::string &attrName, int value){

    // Create a dataset attribute.
    hsize_t dims[1] = { 1 };
    H5::DataSpace attr_dataspace = H5::DataSpace (1, dims);
    H5::Attribute attribute = group->createAttribute( attrName.c_str(), H5::PredType::STD_I32BE,
            attr_dataspace);

    // Write the attribute data.
    int attr_data[1] = {value};
    attribute.write( H5::PredType::NATIVE_INT, attr_data);
}