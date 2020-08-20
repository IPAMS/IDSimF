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
 PSim_interpolatedField.hpp

 Three dimensional field (scalars and 3d vectors) with spatial interpolation of the field values

 ****************************/

#ifndef BTree_interpolatedField_hpp
#define BTree_interpolatedField_hpp

#include <vector>
#include <array>
#include <string>

//forward declare own classes:
namespace Core{
    class Vector;
}

namespace ParticleSimulation{


    /**
     *  Three dimensional rectangular and regular field, capable of containing multiple spatially resolved scalar and
     *  3d vector fields and allowing spatial interpolation of the field values
     */

    class InterpolatedField {
        // TODO: implement an additional equidistant / regular grid mode, with much simplified
        // calculation of the indices of the nodes surrounding an spatial position within the grid

    public:
        explicit InterpolatedField(const std::string &hdf5Filename);

        double getScalar(int ix, int iy, int iz, int fieldIndex);
        double getInterpolatedScalar(double x,double y, double z, int fieldIndex);
        Core::Vector getVector(int ix, int iy, int iz, int fieldIndex);
        Core::Vector getInterpolatedVector(double x,double y, double z, int fieldIndex);
        
        std::vector<std::vector<double>> getGrid();
        std::array<double,6> getBounds();
        std::array<int ,3> findLowerBoundIndices(double x, double y, double z);

    private:
        std::vector<double> gridPointsX_; ///< points of the spatial grid in X direction
        std::vector<double> gridPointsY_; ///< points of the spatial grid in Y direction
        std::vector<double> gridPointsZ_; ///< points of the spatial grid in Z direction
        int gridDimensions_[3];           ///< numbers of spatial grid points in the directions
        std::array<double,6> bounds_;     ///< lower and upper bounds of the spatial grid in the directions

        //std::vector<int> fieldsComponents_;            ///< number of components of the individual data fields
        std::vector<std::string> fieldNames_;
        std::vector<bool> isVector_; ///< Flag if data field is a vector
        std::vector<std::vector<double>> linearizedFields_;

        void updateBounds_();
        inline void enforceScalar_(int fieldIndex);
        int linearizedIndexScalar_(int ix, int iy, int iz);
        int linearizedIndexVector_(int ix, int iy, int iz);
        template<typename datT> datT interpolate_(double x, double y, double z, int fieldIndex);

    };
}

#endif /* BTree_interpolatedFieldEigen_hpp */
