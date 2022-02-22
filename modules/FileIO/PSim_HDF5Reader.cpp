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

#include "PSim_HDF5Reader.hpp"

ParticleSimulation::HDF5Reader::HDF5Reader(const std::string &hdf5Filename) {
    h5f_ = std::make_unique<H5::H5File>(hdf5Filename.c_str(), H5F_ACC_RDONLY);
}

hsize_t ParticleSimulation::HDF5Reader::numberOfObjectsInGroup(std::string groupName) const{
    H5::Group group (h5f_->openGroup(groupName.c_str()));
    return group.getNumObjs();
}


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
    H5O_info_t object_info;
    H5Oget_info(object, &object_info);

    //Write object name to vector if it is a dataset:
    if (object_info.type == H5O_TYPE_DATASET){
        auto *nameVec =  static_cast<std::vector<std::string>*>(opdata);
        nameVec->emplace_back(std::string(name));
    }
    H5Oclose(object);
    return 0;
}

std::vector<std::string> ParticleSimulation::HDF5Reader::namesOfObjectsInGroup(std::string groupName) const{
    H5::Group group = h5f_->openGroup(groupName.c_str());
    std::vector<std::string> objectNames;
    H5Literate(group.getId(), H5_INDEX_NAME, H5_ITER_INC, nullptr, collectObjectNames, &objectNames);
    return objectNames;
}

std::vector<std::string> ParticleSimulation::HDF5Reader::namesOfDatasetsInGroup(std::string groupName) const{
    H5::Group group = h5f_->openGroup(groupName.c_str());
    std::vector<std::string> objectNames;
    H5Literate(group.getId(), H5_INDEX_NAME, H5_ITER_INC, nullptr, collectDatasetNames, &objectNames);
    return objectNames;
}

int ParticleSimulation::HDF5Reader::datasetNDims(std::string datasetName) const{
    H5::DataSet ds = h5f_->openDataSet(datasetName.c_str());
    return ds.getSpace().getSimpleExtentNdims();
}