#include <type_traits>
#include <sstream>
#include <iostream>

template <hsize_t NDIMS, typename DTYPE>
DTYPE ParticleSimulation::HDF5Reader::DataField<NDIMS, DTYPE>::get(std::array<hsize_t, NDIMS> indices){
    hsize_t linIndex;
    if (NDIMS == 4) {
        linIndex = indices[0] * dims[3] * dims[2] * dims[1] + indices[1] * dims[3] * dims[2] + indices[2] * dims[3] + indices[3];
    }
    else if (NDIMS == 3) {
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

template <hsize_t NDIMS>
ParticleSimulation::HDF5Reader::DataField<NDIMS, double>
    ParticleSimulation::HDF5Reader::readDataset(std::string datasetName)
{
    //get the dataspace
    H5::DataSet ds = h5f_->openDataSet(datasetName.c_str());
    return readDataset_<NDIMS>(ds);
}

//template <hsize_t NDIMS>
//ParticleSimulation::HDF5Reader::DataField<NDIMS, double>


/**
 * Reads an one dimensional attribute vector from the HDF5 file
 *
 * @tparam DTYPE The type of the attribute vector
 * @param groupName The name of the group the attribute vector is in
 * @param attributeName The name of the attribute to read
 * @return A std::vector with the data from the attribute vector
 */
template<typename DTYPE>
std::vector<DTYPE> ParticleSimulation::HDF5Reader::readAttributeVector(std::string groupName, std::string attributeName) {

    H5::Group group (h5f_->openGroup(groupName.c_str()));
    H5::Attribute attr(group.openAttribute(attributeName.c_str()));
    H5::DataSpace dataspace = attr.getSpace();

    //get dimensions:
    hsize_t dims[1];
    int nDims = dataspace.getSimpleExtentDims(dims, NULL);
    if (nDims != 1){
        std::stringstream ss;
        ss << "Attribute " << attributeName <<" is not a one dimensional attribute vector";
        throw (std::invalid_argument(ss.str()));
    }

    std::vector<DTYPE> result;

    if constexpr(std::is_same<DTYPE, std::string>::value) {
        char** datBuf = new char* [dims[0]];

        H5::StrType strdatatype(H5::PredType::C_S1, H5T_VARIABLE); // of length 256 characters
        attr.read(strdatatype, datBuf);
        for (int i = 0; i<dims[0]; ++i) {
            result.emplace_back(datBuf[i]);
        }
    } else if constexpr(std::is_same<DTYPE, int>::value || std::is_same<DTYPE, double>::value) {
            const H5::PredType* datType;

            if constexpr(std::is_same<DTYPE, int>::value){
                datType = &H5::PredType::NATIVE_INT;
            } else if constexpr(std::is_same<DTYPE, double>::value){
                datType = &H5::PredType::NATIVE_DOUBLE;
            }

            DTYPE datBuf[dims[0]];
            attr.read(*datType, datBuf);
            for (int i=0; i<dims[0]; ++i){
                result.emplace_back(datBuf[i]);
            }
    } else {
        std::stringstream ss;
        ss << "Reading of " << attributeName <<" with illegal datatype";
        throw (std::invalid_argument(ss.str()));
    }

    return result;
}

template <hsize_t NDIMS>
ParticleSimulation::HDF5Reader::DataField<NDIMS, double>
ParticleSimulation::HDF5Reader::readDataset_(H5::DataSet ds)
{
    //get the dataspace
    H5::DataSpace dataspace = ds.getSpace();

    //get dimensions:
    hsize_t dims[NDIMS];
    int nDims = dataspace.getSimpleExtentDims(dims,NULL);

    //prepare return object and prepare to read from HDF5 file:
    DataField<NDIMS,double> dField;
    dField.rank = NDIMS;
    hsize_t nElements = 1;
    hsize_t offset[NDIMS];

    for (int i=0; i<NDIMS; ++i){
        nElements *= dims[i];
        dField.dims[i] = dims[i];
        offset[i] = 0;
    }

    //define selected hyperslab:
    dataspace.selectHyperslab(H5S_SELECT_SET, dims, offset);

    //define memory dataspace and hyperslab:
    H5::DataSpace memspace(NDIMS, dims);
    memspace.selectHyperslab(H5S_SELECT_SET, dims, offset);

    //read:
    const H5::PredType* nativeType;
    nativeType = &H5::PredType::NATIVE_DOUBLE;

    // init a new stl::vector and write the data directly to it:
    dField.data = std::vector<double>(nElements);
    double* writeBuf = dField.data.data(); //direct access to underlying c style array
    ds.read(writeBuf, *nativeType, memspace, dataspace);
    return dField;
}

