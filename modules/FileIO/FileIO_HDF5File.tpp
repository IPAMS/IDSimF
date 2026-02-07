#include <type_traits>
#include <sstream>
#include <iostream>

template <hsize_t NDIMS, typename DTYPE>
DTYPE FileIO::HDF5File::DataField<NDIMS, DTYPE>::get(std::array<hsize_t, NDIMS> indices){
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

template<typename DTYPE>
void FileIO::HDF5File::writeDatasetAttribute(H5::DataSet &dataSet, const std::string& attrName, const std::vector<DTYPE> &values) {
    hsize_t nVals = values.size();
    hsize_t dims[1] = { nVals };
    H5::DataSpace attrDataSpace = H5::DataSpace (1, dims);

    if constexpr(std::is_same_v<DTYPE, std::string>) {
        H5::StrType strdatatype(H5::PredType::C_S1, H5T_VARIABLE); // of length 256 characters

        //Due to a abi bug between the std lib and hdf, we need to provide raw strings to the hdf5 methods
        //Save array of pointers to the raw c strings and use that array as data buffer
        std::vector<const char*> dataBuf(nVals);
        const char** dataBufArray = dataBuf.data();
        for (std::size_t i = 0; i < nVals; ++i)
        {
            dataBufArray[i] = values[i].c_str();
        }
        H5::Attribute stringAttribute = dataSet.createAttribute(attrName.c_str(), strdatatype, attrDataSpace);
        stringAttribute.write(strdatatype, dataBufArray);
    }
    else if constexpr(std::is_same_v<DTYPE, int>) {
        H5::Attribute intAttribute = dataSet.createAttribute(attrName.c_str(), H5::PredType::STD_I32BE, attrDataSpace);
        intAttribute.write(H5::PredType::NATIVE_INT, values.data());
    }
    else if constexpr(std::is_same_v<DTYPE, double>) {
        H5::Attribute doubleAttribute = dataSet.createAttribute(attrName.c_str(), H5::PredType::IEEE_F32BE, attrDataSpace);
        doubleAttribute.write(H5::PredType::NATIVE_DOUBLE, values.data());
    }
    else {
        //use workaround since static_assert(false) leads to compiler / template instantiation problems
        static_assert(!std::is_same_v<DTYPE, DTYPE>);
    }
}

template <hsize_t NDIMS>
FileIO::HDF5File::DataField<NDIMS, double>
    FileIO::HDF5File::readDataset(std::string datasetName) const
{
    //get the dataspace
    H5::DataSet ds = h5f_->openDataSet(datasetName.c_str());
    return readDataset_<NDIMS>(ds);
}

/**
 * Reads an one dimensional attribute vector from the HDF5 file
 *
 * @tparam DTYPE The type of the attribute vector
 * @param groupName The name of the group the attribute vector is in
 * @param attributeName The name of the attribute to read
 * @return A std::vector with the data from the attribute vector
 */
template<typename DTYPE>
std::vector<DTYPE> FileIO::HDF5File::readAttributeVector(std::string groupName, std::string attributeName) const{

    H5::Group group (h5f_->openGroup(groupName.c_str()));
    H5::Attribute attr(group.openAttribute(attributeName.c_str()));
    H5::DataSpace dataspace = attr.getSpace();

    //get dimensions:
    hsize_t dims[1];
    int nDims = dataspace.getSimpleExtentDims(dims, nullptr);
    if (nDims != 1){
        std::stringstream ss;
        ss << "Attribute " << attributeName <<" is not a one dimensional attribute vector";
        throw (std::invalid_argument(ss.str()));
    }

    std::vector<DTYPE> result;

    if constexpr(std::is_same_v<DTYPE, std::string>) {
        char** datBuf = new char* [dims[0]];

        H5::StrType strdatatype(H5::PredType::C_S1, H5T_VARIABLE);
        attr.read(strdatatype, datBuf);
        for (hsize_t i = 0; i<dims[0]; ++i) {
            result.emplace_back(datBuf[i]);
        }
    } else if constexpr(std::is_same_v<DTYPE, int> || std::is_same_v<DTYPE, double>) {
            const H5::PredType* datType;

            if constexpr(std::is_same_v<DTYPE, int>){
                datType = &H5::PredType::NATIVE_INT;
            } else if constexpr(std::is_same_v<DTYPE, double>){
                datType = &H5::PredType::NATIVE_DOUBLE;
            }

            std::vector<DTYPE> datBuf(dims[0]);
            attr.read(*datType, datBuf.data());
            for (hsize_t i=0; i<dims[0]; ++i){
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
FileIO::HDF5File::DataField<NDIMS, double>
FileIO::HDF5File::readDataset_(H5::DataSet ds) const
{
    //get the dataspace
    H5::DataSpace dataspace = ds.getSpace();

    //get dimensions:
    hsize_t dims[NDIMS];
    dataspace.getSimpleExtentDims(dims, nullptr);

    //prepare return object and prepare to read from HDF5 file:
    DataField<NDIMS,double> dField;
    dField.rank = NDIMS;
    hsize_t nElements = 1;
    hsize_t offset[NDIMS];

    for (hsize_t i=0; i<NDIMS; ++i){
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

