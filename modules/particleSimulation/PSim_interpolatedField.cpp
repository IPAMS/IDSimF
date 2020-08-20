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

#include "PSim_interpolatedField.hpp"
#include "PSim_HDF5Reader.hpp"
#include "Core_vector.hpp"
#include <algorithm>
#include <iterator>
#include <exception>

//TODO: add comments

template <hsize_t dims> using dataField =
    ParticleSimulation::HDF5Reader::DataField<dims, double>;


ParticleSimulation::InterpolatedField::InterpolatedField(const std::string &hdf5Filename) {
    ParticleSimulation::HDF5Reader h5Reader(hdf5Filename);

    dataField<1> gridPointsXDf =  h5Reader.readDataset<1>("/grid_points/x");
    gridPointsX_ = gridPointsXDf.data;

    dataField<1> gridPointsYDf =  h5Reader.readDataset<1>("/grid_points/y");
    gridPointsY_ = gridPointsYDf.data;

    dataField<1> gridPointsZDf =  h5Reader.readDataset<1>("/grid_points/z");
    gridPointsZ_ = gridPointsZDf.data;

    updateBounds_();

    gridDimensions_[0] = gridPointsX_.size();
    gridDimensions_[1] = gridPointsY_.size();
    gridDimensions_[2] = gridPointsZ_.size();

    // retrieve names of fields in HDF5 file
    fieldNames_ = h5Reader.namesOfDatasetsInGroup("/fields");

    //For loop to iterate through all values in a vector / array:
    for(const auto &fieldName: fieldNames_) {
        //std::cout << temp <<std::endl;
        //h5Reader.readDataset<>()
        std::string fullFieldName = "/fields/"+fieldName;
        const int nDims = h5Reader.datasetNDims(fullFieldName);
        if (nDims==3) {
            //scalar field
            isVector_.push_back(false);
        }
        else if (nDims==4) {
            isVector_.push_back(true);
        }
        else {
            std::stringstream ss;
            ss << "Dataset " << fieldName << " has illegal dimensionality";
            throw (std::invalid_argument(ss.str()));
        }

        linearizedFields_.emplace_back(std::vector<double>());
        unsigned int counter = 0;

        // Read field data
        if (isVector_.back()) {
            dataField<4> df = h5Reader.readDataset<4>(fullFieldName);

            //TODO: throw exception if 4th dimension (the vector components dimension) is not 3 (not a 3d vector)


            std::array<hsize_t, 4> dfIndices = {0, 0, 0, 0};
            std::vector<double>& lf = linearizedFields_.back();

            for (int zi = 0; zi<gridDimensions_[2]; ++zi) {
                dfIndices[2] = zi;
                for (int yi = 0; yi<gridDimensions_[1]; ++yi) {
                    dfIndices[1] = yi;
                    for (int xi = 0; xi<gridDimensions_[0]; ++xi) {
                        dfIndices[0] = xi;

                        //4th dimension is the vector components dimension:
                        dfIndices[3] = 0;
                        lf.emplace_back(df.get(dfIndices));
                        dfIndices[3] = 1;
                        lf.emplace_back(df.get(dfIndices));
                        dfIndices[3] = 2;
                        lf.emplace_back(df.get(dfIndices));

                        ++counter;
                    }
                }
            }

        }
        else {
            dataField<3> df = h5Reader.readDataset<3>(fullFieldName);
            std::array<hsize_t, 3> dfIndices = {0, 0, 0};
            std::vector<double>& lf = linearizedFields_.back();

            for (int zi = 0; zi<gridDimensions_[2]; ++zi) {
                dfIndices[2] = zi;
                for (int yi = 0; yi<gridDimensions_[1]; ++yi) {
                    dfIndices[1] = yi;
                    for (int xi = 0; xi<gridDimensions_[0]; ++xi) {
                        dfIndices[0] = xi;
                        lf.emplace_back(df.get(dfIndices));
                        ++counter;
                    }
                }
            }
        }
    }
}

double ParticleSimulation::InterpolatedField::getScalar(int ix, int iy, int iz, int fieldIndex) {
    enforceScalar_(fieldIndex);
    return linearizedFields_[fieldIndex][linearizedIndexScalar_(ix,iy,iz)];
}

double ParticleSimulation::InterpolatedField::getInterpolatedScalar(double x, double y, double z, int fieldIndex) {
    enforceScalar_(fieldIndex);
    return interpolate_<double>(x, y, z, fieldIndex);
}

Core::Vector ParticleSimulation::InterpolatedField::getVector(int ix, int iy, int iz, int fieldIndex) {
    if (!isVector_[fieldIndex]){
        std::stringstream ss;
        ss << "Data field " << fieldIndex <<" is not a vector field";
        throw (std::invalid_argument(ss.str()));
    }

    int vectorIndex = linearizedIndexVector_(ix, iy, iz);
    return {
        linearizedFields_[fieldIndex][vectorIndex],
        linearizedFields_[fieldIndex][vectorIndex+1],
        linearizedFields_[fieldIndex][vectorIndex+2]
    };
}

Core::Vector ParticleSimulation::InterpolatedField::getInterpolatedVector(double x, double y, double z, int fieldIndex) {
    if (!isVector_[fieldIndex]){
        std::stringstream ss;
        ss << "Data field " << fieldIndex <<" is not a vector field";
        throw (std::invalid_argument(ss.str()));
    }
    auto result = interpolate_<std::array<double, 3>>(x, y, z, fieldIndex);
    return Core::Vector(result[0], result[1], result[2]);
}

std::vector<std::vector<double>> ParticleSimulation::InterpolatedField::getGrid() {
    std::vector<std::vector<double>> result;
    result.push_back(gridPointsX_);
    result.push_back(gridPointsY_);
    result.push_back(gridPointsZ_);
    return result;
}

std::array<double,6> ParticleSimulation::InterpolatedField::getBounds() {
    return bounds_;
}

std::array<int,3> ParticleSimulation::InterpolatedField::findLowerBoundIndices(double x, double y, double z) {
    auto lowerX = std::lower_bound(gridPointsX_.begin(), gridPointsX_.end(), x);
    int indexLowerX = lowerX-gridPointsX_.begin();
    auto lowerY = std::lower_bound(gridPointsY_.begin(), gridPointsY_.end(), y);
    int indexLowerY = lowerY-gridPointsY_.begin();
    auto lowerZ = std::lower_bound(gridPointsZ_.begin(), gridPointsZ_.end(), z);
    int indexLowerZ = lowerZ-gridPointsZ_.begin();

    std::array<int,3> indexes = {indexLowerX,indexLowerY,indexLowerZ};

    return indexes;
}

void ParticleSimulation::InterpolatedField::updateBounds_() {
    bounds_[0] = gridPointsX_.front();
    bounds_[1] = gridPointsX_.back();
    bounds_[2] = gridPointsY_.front();
    bounds_[3] = gridPointsY_.back();
    bounds_[4] = gridPointsZ_.front();
    bounds_[5] = gridPointsZ_.back();
}

void ParticleSimulation::InterpolatedField::enforceScalar_(int fieldIndex) {
    if (isVector_[fieldIndex]){
        std::stringstream ss;
        ss << "Data field " << fieldIndex <<" is not a scalar field";
        throw (std::invalid_argument(ss.str()));
    }
}

int ParticleSimulation::InterpolatedField::linearizedIndexScalar_(int ix, int iy, int iz) {
    int index = iz*gridDimensions_[1]*gridDimensions_[0] +
                iy*gridDimensions_[0] +
                ix;
    return index;
}

int ParticleSimulation::InterpolatedField::linearizedIndexVector_(int ix, int iy, int iz) {
    int index = iz * gridDimensions_[1] * gridDimensions_[0] * 3 +
                iy * gridDimensions_[0] * 3 +
                ix * 3;

    return index;
}

template<typename datT>
datT ParticleSimulation::InterpolatedField::interpolate_(double x, double y, double z, int fieldIndex) {

    //if out of bounds: just return 0 (as the VTK interpolator also does)
    //FIXME: throw exception here
    if (x<= bounds_[0] || x>= bounds_[1] || y<= bounds_[2] || y>= bounds_[3] || z<= bounds_[4] || z>= bounds_[5]){
        std::stringstream ss;
        ss << "Point with coordinates " << x <<" "<< y <<" "<< z <<" is not in bounds of interpolated field";
        throw (std::invalid_argument(ss.str()));
    }

    std::vector<double>* field = &(linearizedFields_[fieldIndex]);

    std::array<int,3> lowerBoundIndices = findLowerBoundIndices(x, y, z);
    std::size_t xiLower = lowerBoundIndices[0]-1;
    int xiUpper = lowerBoundIndices[0];
    int yiLower = lowerBoundIndices[1]-1;
    int yiUpper = lowerBoundIndices[1];
    int ziLower = lowerBoundIndices[2]-1;
    int ziUpper = lowerBoundIndices[2];

    double xLower = gridPointsX_.at(xiLower);
    double xUpper = gridPointsX_.at(xiUpper);
    double yLower = gridPointsY_.at(yiLower);
    double yUpper = gridPointsY_.at(yiUpper);
    double zLower = gridPointsZ_.at(ziLower);
    double zUpper = gridPointsZ_.at(ziUpper);

    double xd = (x-xLower) /(xUpper - xLower);
    double yd = (y-yLower) /(yUpper - yLower);
    double zd = (z-zLower) /(zUpper - zLower);

    if constexpr (std::is_same_v<datT, double>) {

        double c_00 = field->operator[](linearizedIndexScalar_(xiLower, yiLower, ziLower))*(1-xd)
                    + field->operator[](linearizedIndexScalar_(xiUpper, yiLower, ziLower))*xd;

        double c_01 = field->operator[](linearizedIndexScalar_(xiLower, yiLower, ziUpper))*(1-xd)
                    + field->operator[](linearizedIndexScalar_(xiUpper, yiLower, ziUpper))*xd;

        double c_10 = field->operator[](linearizedIndexScalar_(xiLower, yiUpper, ziLower))*(1-xd)
                    + field->operator[](linearizedIndexScalar_(xiUpper, yiUpper, ziLower))*xd;

        double c_11 = field->operator[](linearizedIndexScalar_(xiLower, yiUpper, ziUpper))*(1-xd)
                    + field->operator[](linearizedIndexScalar_(xiUpper, yiUpper, ziUpper))*xd;

        double c_0 = c_00*(1-yd) + c_10*yd;
        double c_1 = c_01*(1-yd) + c_11*yd;

        double result = c_0*(1-zd) + c_1*zd;

        return result;
    }
    else if constexpr (std::is_same_v<datT, std::array<double, 3>>) {

        int i_000 = linearizedIndexVector_(xiLower, yiLower, ziLower);
        int i_001 = linearizedIndexVector_(xiLower, yiLower, ziUpper);
        int i_010 = linearizedIndexVector_(xiLower, yiUpper, ziLower);
        int i_011 = linearizedIndexVector_(xiLower, yiUpper, ziUpper);
        int i_100 = linearizedIndexVector_(xiUpper, yiLower, ziLower);
        int i_101 = linearizedIndexVector_(xiUpper, yiLower, ziUpper);
        int i_110 = linearizedIndexVector_(xiUpper, yiUpper, ziLower);
        int i_111 = linearizedIndexVector_(xiUpper, yiUpper, ziUpper);

        std::array<double, 3> result{0.0,0.0,0.0};

        for (int i=0; i<3; ++i){
            double c_00 = field->operator[](i_000+i)*(1-xd) + field->operator[](i_100+i)*xd;
            double c_01 = field->operator[](i_001+i)*(1-xd) + field->operator[](i_101+i)*xd;
            double c_10 = field->operator[](i_010+i)*(1-xd) + field->operator[](i_110+i)*xd;
            double c_11 = field->operator[](i_011+i)*(1-xd) + field->operator[](i_111+i)*xd;

            double c_0 = c_00*(1-yd) + c_10*yd;
            double c_1 = c_01*(1-yd) + c_11*yd;

            result[i] = c_0*(1-zd) + c_1*zd;

        }

        return result;
    }
}