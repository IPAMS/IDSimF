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
 PSim_averageChargePositionWriter.hpp

 File writer for average position of the net charge in a Barnes Hut Tree

 ****************************/

#ifndef BTree_average_charge_position_writer_hpp
#define BTree_average_charge_position_writer_hpp

#include <fstream>
#include <string>

//forward declare own classes:
namespace BTree{
    class Tree;
}

namespace ParticleSimulation {
    /**
     * Filewriter to write the average position of the net charge in a Core to a file
     */
    class AverageChargePositionWriter{
    public:
        explicit AverageChargePositionWriter(std::string transientFilename);
        ~AverageChargePositionWriter();
        void writeTimestep(const BTree::Tree& tree, double time);

    private:
        std::ofstream* transientFile_; ///< file handle of the file to write to

    };
}

#endif /* BTree_average_charge_position_writer_hpp */
