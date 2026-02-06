/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2026 - Physical and Theoretical Chemistry /
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
 test_hdf5_util.hpp

 Description

 ****************************/
#ifndef IDSIMF_TEST_HDF5_UTIL_HPP
#define IDSIMF_TEST_HDF5_UTIL_HPP

#include "H5Cpp.h"

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
    CHECK(NDIMS == rank);

    //get dimensions:
    hsize_t dims[NDIMS];
    hsize_t nDims = (hsize_t)dataspace.getSimpleExtentDims(dims, nullptr);
    CHECK(NDIMS == nDims);

    //prepare return object and prepare to read from HDF5 file:
    DataField<NDIMS,double> dField;
    dField.rank = nDims;
    hsize_t nElements = 1;
    std::vector<hsize_t> offset(nDims);
    //hsize_t offset[nDims];
    //hsize_t count[nDims];

    for (hsize_t i=0; i<NDIMS; ++i){
        nElements *= dims[i];
        dField.dims[i] = dims[i];
        offset[i] = 0;
        //count[i] = dims[i];
    }

    //define selected hyperslab:
    dataspace.selectHyperslab(H5S_SELECT_SET, dims, offset.data());

    //define memory dataspace and hyperslab:
    H5::DataSpace memspace(static_cast<int>(NDIMS), dims);
    memspace.selectHyperslab(H5S_SELECT_SET, dims, offset.data());

    //read:
    std::vector<double> datBuf(nElements);
    double* datBufArray = datBuf.data();
    const H5::PredType* nativeType;
    nativeType = &H5::PredType::NATIVE_DOUBLE;
    ds.read(datBufArray,*nativeType,memspace,dataspace);

    for (hsize_t i=0; i<nElements; ++i){
        dField.data.emplace_back(datBufArray[i]);
    }

    return dField;
}

inline std::vector<std::string> readStringAttribute(H5::Group& group, std::string attrName){
    H5::Attribute attr(group.openAttribute(attrName.c_str()));
    H5::DataSpace dataspace = attr.getSpace();

    //get dimensions:
    hsize_t dims[1];
    int nDims = dataspace.getSimpleExtentDims(dims, nullptr);
    CHECK(nDims == 1);
    std::vector<std::string> result;
    char **datBuf = new char*[dims[0]];
    H5::StrType strdatatype(H5::PredType::C_S1, H5T_VARIABLE); // of length 256 characters
    attr.read(strdatatype, datBuf);
    for (hsize_t i = 0; i < dims[0]; ++i) {
        result.emplace_back(datBuf[i]);
    }
    return (result);
}

// Todo: Compact code for integral types with templates
inline std::vector<int> readIntAttribute(H5::Group& group, std::string attrName){
    H5::Attribute attr(group.openAttribute(attrName.c_str()));
    H5::DataSpace dataspace = attr.getSpace();

    //get dimensions:
    hsize_t dims[1];
    int nDims = dataspace.getSimpleExtentDims(dims, nullptr);
    CHECK(nDims == 1);
    std::vector<int> result(dims[0]);
    attr.read(H5::PredType::NATIVE_INT, result.data());
    return result;
}

inline std::vector<int> readHSizetAttribute(H5::Group& group, std::string attrName){
    H5::Attribute attr(group.openAttribute(attrName.c_str()));
    H5::DataSpace dataspace = attr.getSpace();

    //get dimensions:
    hsize_t dims[1];
    int nDims = dataspace.getSimpleExtentDims(dims, nullptr);
    CHECK(nDims == 1);
    std::vector<int> result(dims[0]);
    attr.read(H5::PredType::NATIVE_UINT64, result.data());
    return result;
}

inline std::vector<double> readDoubleAttribute(H5::Group& group, std::string attrName){
    H5::Attribute attr(group.openAttribute(attrName.c_str()));
    H5::DataSpace dataspace = attr.getSpace();

    //get dimensions:
    hsize_t dims[1];
    int nDims = dataspace.getSimpleExtentDims(dims, nullptr);
    CHECK(nDims == 1);
    std::vector<double> result(dims[0]);
    attr.read(H5::PredType::NATIVE_DOUBLE, result.data());
    return result;
}

#endif //IDSIMF_TEST_HDF5_UTIL_HPP