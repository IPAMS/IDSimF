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
 PSim_ionCloudReader.hpp

 Reader class for structured simple ion clouds / ion initialization data from csv files

 ****************************/

#ifndef PSim_ionCloudReader_hpp
#define PSim_ionCloudReader_hpp

#include <stdexcept>
#include <string>
#include <vector>
#include <memory>

//forward declare own classes:
namespace BTree{
    class Particle;
}

namespace ParticleSimulation{

    /**
     * Individual exception class for problems in ion cloud files
     */
    class IonCloudFileException : public std::runtime_error {
        public:
            explicit IonCloudFileException (const std::string msg): std::runtime_error(msg) {}
    };

    /**
     * File reader for importing structured simple ion clouds / ion initialization data from files
     */
    class IonCloudReader {
        public:
            std::vector<std::unique_ptr<BTree::Particle>> readIonCloud(std::string filename);

        private:
            const std::string delimiter_ = ";"; ///<A delimiter for the columns in the ion cloud files
    };
}

#endif /* BTree_vtkFieldReader_hpp */
