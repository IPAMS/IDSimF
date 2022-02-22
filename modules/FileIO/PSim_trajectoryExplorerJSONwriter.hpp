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

 Writer class for JSON files containing ion trajectories

 ****************************/

#ifndef BTree_trajectoryExplorerJSONwriter_hpp
#define BTree_trajectoryExplorerJSONwriter_hpp

#include "PSim_trajectoryWriterDefs.hpp"
#include <ostream>
#include <vector>
#include <functional>
#include <memory>

//forward declare own classes:
namespace BTree{
    class Tree;
    class Particle;
}

namespace FileIO {

    /**
     * File writer class for JSON files containing ion trajectories.
     *
     * The spatial dimensions and the time dimension of the data to write can be scaled.
     * (This writer assumes, that the number of ions in the simulation does not change)
     */
    class TrajectoryExplorerJSONwriter {

    public:
        explicit TrajectoryExplorerJSONwriter(std::string jsonFilename);
        ~TrajectoryExplorerJSONwriter();

        void setScales(double scale, double timeScale);

        void writeTimestep(std::vector<BTree::Particle*>& particles, double time, bool lastTime);

        void writeTimestep(
                std::vector<BTree::Particle*>& particles,
                const partAttribTransformFctType& particleParameterTransformFct,
                double time, bool lastTime);

        void writeTimestep(
                std::vector<BTree::Particle*>& particles,
                const partAttribTransformFctType& particleParameterTransformFct,
                std::vector<additionalParamPair> timestepAdditionalParameters,
                double time, bool lastTime);

        void writeSplatTimes(std::vector<BTree::Particle*>& particles);
        void writeIonMasses(std::vector<BTree::Particle*>& particles);

    private:
        std::unique_ptr<std::ofstream> jsonFile_; ///< file stream to write to
        double scale_ = 1.0; ///< spatial scaling factor
        double timeScale_ = 1.0; ///< time scaling factor

        void initFile_();
        void closeFile_();

        void writeIonPosition_(BTree::Particle* particle);
        void writeIonSplatTime_(BTree::Particle* particle, bool lastParticle);
    };
}

#endif /* BTree_trajectoryExplorerJSONwriter_hpp */
