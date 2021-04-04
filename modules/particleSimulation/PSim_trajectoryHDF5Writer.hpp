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
 PSim_trajectoryExplorerJSONwriter.hpp

 Writer class for HDF5 files containing ion trajectories

 ****************************/

#ifndef BTree_trajectoryHDF5Writer_hpp
#define BTree_trajectoryHDF5Writer_hpp

#include "BTree_particle.hpp"
#include "PSim_trajectoryWriterDefs.hpp"
#include "H5Cpp.h"
#include <string>
#include <vector>
#include <memory>


namespace BTree{
    class Particle;
}

namespace ParticleSimulation{

    class ParticleStartSplatTracker;

    class TrajectoryHDF5Writer{
    public:
        TrajectoryHDF5Writer(const std::string &hdf5Filename, bool compression = true);

        TrajectoryHDF5Writer(const std::string &hdf5Filename,
                             std::vector<std::string> auxParamsNames,
                             additionalPartParamFctType particleParameterTransformFct, bool compression = true);

        void writeTimestep(std::vector<BTree::Particle*>& particles, double time);

        void writeVectorDataset(std::string dsName, std::vector<double> values);
        void writeTrajectoryAttribute(std::string attrName, int value);
        void writeTrajectoryAttribute(std::string attrName, std::vector<double> values);
        void writeTrajectoryAttribute(std::string attrName, std::vector<std::string> values);
        void finalizeTrajectory();
        void writeSplatTimes(std::vector<BTree::Particle*> &particles);
        void writeStartSplatData(ParticleStartSplatTracker tracker);

    private:
        constexpr int static FILE_TYPE_VERSION = 3;   ///<File type version identifier of the files written

        void writeAuxTimestep_(std::vector<BTree::Particle*> &particles);
        void writeAttribute_(std::unique_ptr<H5::Group>& group, const std::string& attrName, int value);

        //int nTimestepsWritten_;
        bool compression_;
        bool isAuxWritten_;

        std::unique_ptr<H5::H5File > h5f_;
        std::unique_ptr<H5::Group> baseGroup_;
        std::unique_ptr<H5::Group> timeStepGroup_;
        std::unique_ptr<H5::DataSet> dsetTimesteps_;

        hsize_t nAuxParams_;          ///< number auf auxiliary parameters for the parameters
        hsize_t sizeTimesteps_[1];    ///< ... and the timesteps
        hsize_t offsetScalarLike_[3]; ///< current offset in the vector of scalars datasets
        hsize_t slabDimsTimestep_[1]; ///< size of one timestep frame in the timesteps

        H5::DataSpace memspaceTimestep_;

        additionalPartParamFctType particleParameterTransformFct_;
    };
}


#endif //BTree_trajectoryHDF5Writer_hpp
