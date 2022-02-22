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

#include "PSim_scalar_writer.hpp"
#include "BTree_tree.hpp"
#include "BTree_particle.hpp"

/**
 * Constructor: Creates a file writer for a given result file
 * @param transientFilename the name of the file to write to
 */
FileIO::Scalar_writer::Scalar_writer(std::string transientFilename){
    transientFile_ = new std::ofstream();
    transientFile_->open(transientFilename);
}

/**
 * Destructor
 */
FileIO::Scalar_writer::~Scalar_writer(){
    if (transientFile_ != nullptr){
        transientFile_->flush();
        delete (transientFile_);
    }
}

/**
 * Writes a scalar value to the result file
 * @param intValue the scalar value to write
 * @param time the time of the current time step
 */
void FileIO::Scalar_writer::writeTimestep(int intValue, double time){
    *transientFile_<<time<<" "<<intValue<<std::endl;
}

/**
 * Writes a scalar value to the result file
 * @param sizeValue the scalar value to write
 * @param time the time of the current time step
 */
void FileIO::Scalar_writer::writeTimestep(std::size_t sizeValue, double time){
    *transientFile_<<time<<" "<<sizeValue<<std::endl;
}

/**
 * Writes a unsigned integer scalar value to the result file
 * @param unsignedIntValue the scalar value to write
 * @param time the time of the current time step
 */
void FileIO::Scalar_writer::writeTimestep(unsigned int unsignedIntValue, double time){
    *transientFile_<<time<<" "<<unsignedIntValue<<std::endl;
}

/**
 * Writes a scalar value to the result file
 * @param doubleValue the scalar value to write
 * @param time the time of the current time step
 */
void FileIO::Scalar_writer::writeTimestep(double doubleValue, double time){
    *transientFile_<<time<<" "<<doubleValue<<std::endl;
}

/**
 * Writes multiple values simultaneously to the result file
 * @param doubleValues vector of values to write
 * @param time the time of the current time step
 */
void FileIO::Scalar_writer::writeTimestep(std::vector<double> doubleValues, double time){
    *transientFile_<<time<<" ";
    for(const auto &val: doubleValues){
        *transientFile_<<val<<" ";
    }
    *transientFile_<<std::endl;
}
