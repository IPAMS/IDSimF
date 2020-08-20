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
 PSim_HDF5Reader.hpp

 Description

 ****************************/

#ifndef IDSIMF_PSIM_HDF5READER_HPP
#define IDSIMF_PSIM_HDF5READER_HPP

#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include "H5Cpp.h"

namespace ParticleSimulation{
    class HDF5Reader {

    public:
        // Define some helping data structures and functions .........
        template <hsize_t NDIMS, typename DTYPE> struct DataField{
            hsize_t rank = 0;
            std::array<hsize_t, NDIMS> dims;
            std::vector<DTYPE> data;
            DTYPE get(std::array<hsize_t, NDIMS> indices);
        };

        explicit HDF5Reader(const std::string &hdf5Filename, bool compression = true);

        template <hsize_t NDIMS> DataField<NDIMS, double>
        readDataset(std::string datasetName);

        template<typename DTYPE>
        std::vector<DTYPE> readAttributeVector(std::string groupName, std::string attributeName);

        hsize_t numberOfObjectsInGroup(std::string groupName);
        std::vector<std::string> namesOfObjectsInGroup(std::string groupName);
        std::vector<std::string> namesOfDatasetsInGroup(std::string groupName);
        int datasetNDims(std::string datasetName);


    private:
        std::unique_ptr<H5::H5File > h5f_;

        template <hsize_t NDIMS> DataField<NDIMS, double>
        readDataset_(H5::DataSet ds);
    };
}

#include "PSim_HDF5Reader.tpp"

#endif //IDSIMF_PSIM_HDF5READER_HPP
