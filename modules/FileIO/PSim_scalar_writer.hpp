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
 PSim_scalar_writer.hpp

 Simple file writer for scalar values

 ****************************/

#ifndef Particle_simulation_scalar_writer
#define Particle_simulation_scalar_writer

#include <fstream>
#include <vector>


//forward declare own classes:
namespace BTree{
    class Vector;
    class Particle;
}

namespace FileIO {
    /**
     * File writer to write scalar values to a result file
     */
    class Scalar_writer{
    public:
        explicit Scalar_writer(std::string transientFilename);
        ~Scalar_writer();

        void writeTimestep(int intValue, double time);
        void writeTimestep(std::size_t sizeValue, double time);
        void writeTimestep(unsigned int unsignedIntValue, double time);
        void writeTimestep(double doubleValue, double time);
        void writeTimestep(std::vector<double> doubleValues, double time);

    private:
        std::ofstream* transientFile_; ///< File stream to the file to write to
    };
}


#endif //Particle_simulation_scalar_writer
