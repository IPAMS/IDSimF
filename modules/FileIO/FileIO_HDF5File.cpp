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

#include "FileIO_HDF5File.hpp"
#include <string>
#include <sstream>
#include <cassert>

FileIO::HDF5File::HDF5File(const std::string &hdf5Filename, FileMode mode){
    if (mode == READ_ONLY) {
        h5f_ = std::make_unique<H5::H5File>(hdf5Filename.c_str(), H5F_ACC_RDONLY);
    }
    else if (mode == WRITE_ONLY) {
        h5f_ = std::make_unique<H5::H5File>(hdf5Filename.c_str(), H5F_ACC_TRUNC);
    }
}

hsize_t FileIO::HDF5File::numberOfObjectsInGroup(std::string groupName) const{
    H5::Group group (h5f_->openGroup(groupName.c_str()));
    return group.getNumObjs();
}


bool FileIO::HDF5File::groupPathExists(std::string groupName) const{
    //traverse full group path and check for existence:
    std::stringstream ss(groupName);
    std::vector<std::string> parts;
    std::string pathPart;

    while (std::getline(ss, pathPart, '/')) {
        if (!pathPart.empty()) {
            parts.push_back(pathPart);
        }
    }

    std::string currentPath = "";
    for (const auto& part : parts) {
        currentPath += "/" + part;
        if (H5Lexists(h5f_->getId(), currentPath.c_str(), H5P_DEFAULT) <=0) {
            return false;
        };
    }
    return true;
}

H5::Group FileIO::HDF5File::createGroup(std::string groupName) const {
    if (!groupPathExists(groupName)) {
        H5::LinkCreatPropList propList;
        propList.setCreateIntermediateGroup(true);
        H5::Group group = h5f_->createGroup(groupName.c_str(), propList);
        return group;
    }
    else {
        H5::Group group = h5f_->openGroup(groupName.c_str());
        return group;
    }
}

H5::DataSet FileIO::HDF5File::initTableDataset(std::string groupName, std::string datasetName, std::size_t nColumns) const{

    //prepare location dataset structures:
    hsize_t dimsDS[2] = {0, nColumns};            // dataset dimensions at creation
    hsize_t maxdimsDS[2] = {H5S_UNLIMITED, nColumns};         // maximum dataset dimensions
    hsize_t chunkDimsDS[2] = {1, nColumns};

    H5::DataSpace dataspaceLocation(2, dimsDS, maxdimsDS);

    H5::DSetCreatPropList propTimesteps;
    propTimesteps.setChunk(2, chunkDimsDS);

    //create actual dataset for times:
    H5::Group group = createGroup(groupName);
    H5::DataSpace dataspaceDS(2,dimsDS,maxdimsDS);
    H5::DataSet dataSet = group.createDataSet(datasetName, H5::PredType::IEEE_F32BE, dataspaceDS, propTimesteps);

    return dataSet;
}

void FileIO::HDF5File::writeRowToTableDataset(H5::DataSet dataSet, const std::vector<double>& data) {

    // extend dataset:
    hsize_t dims[2];
    dataSet.getSpace().getSimpleExtentDims(dims);

    if (dims[1] != data.size()) {
        std::stringstream ss;
        ss << "Invalid input vector for writing a row in HDF5 table. Table columns: "<<dims[1]<< " input vector length: "<<data.size() <<std::endl;
        throw std::invalid_argument(ss.str());
    }

    hsize_t offset[2] = {dims[0], 0};
    dims[0] += 1;
    dataSet.extend(dims);
    H5::DataSpace fileSpace = dataSet.getSpace();

    //write to dataset:
    hsize_t memSpaceDims[2] = {1,dims[1]};
    fileSpace.selectHyperslab(H5S_SELECT_SET, memSpaceDims, offset);
    H5::DataSpace memspace(2, memSpaceDims);
    dataSet.write(data.data(), H5::PredType::NATIVE_DOUBLE, memspace, fileSpace);
}

// some functions to be used on data set iterations:
herr_t collectObjectNames(hid_t /*loc_id*/, const char *name, const H5L_info_t* /*linfo*/, void *opdata)
{
    auto *nameVec =  static_cast<std::vector<std::string>*>(opdata);
    nameVec->emplace_back(std::string(name));

    return 0;
}

herr_t collectDatasetNames(hid_t loc_id, const char *name, const H5L_info_t* /*linfo*/, void *opdata)
{
    // Open the object using its name.
    hid_t object = H5Oopen(loc_id, name, H5P_DEFAULT);

    #ifdef OLD_HDF5_API
        H5O_info_t object_info;
        H5Oget_info(object, &object_info);
    #else
        H5O_info1_t object_info;
        H5Oget_info1(object, &object_info);
    #endif

    //Write object name to vector if it is a dataset:
    if (object_info.type == H5O_TYPE_DATASET){
        auto *nameVec =  static_cast<std::vector<std::string>*>(opdata);
        nameVec->emplace_back(std::string(name));
    }
    H5Oclose(object);
    return 0;
}
// ----------------------------------

std::vector<std::string> FileIO::HDF5File::namesOfObjectsInGroup(std::string groupName) const{
    H5::Group group = h5f_->openGroup(groupName.c_str());
    std::vector<std::string> objectNames;
    H5Literate(group.getId(), H5_INDEX_NAME, H5_ITER_INC, nullptr, collectObjectNames, &objectNames);
    return objectNames;
}

std::vector<std::string> FileIO::HDF5File::namesOfDatasetsInGroup(std::string groupName) const{
    H5::Group group = h5f_->openGroup(groupName.c_str());
    std::vector<std::string> objectNames;
    H5Literate(group.getId(), H5_INDEX_NAME, H5_ITER_INC, nullptr, collectDatasetNames, &objectNames);
    return objectNames;
}

int FileIO::HDF5File::datasetNDims(std::string datasetName) const{
    H5::DataSet ds = h5f_->openDataSet(datasetName.c_str());
    return ds.getSpace().getSimpleExtentNdims();
}