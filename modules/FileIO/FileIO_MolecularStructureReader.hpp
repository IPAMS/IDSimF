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
 FileIO_molecularStructureReader.hpp

 Reader class for molecular structure data from csv files

 ****************************/

#ifndef IDSIM_FILEIO_MOLECULARSTRUCTUREREADER
#define IDSIM_FILEIO_MOLECULARSTRUCTUREREADER

#include <string>

namespace FileIO
{
    class MolecularStructureReader{

        public:
            void readMolecularStructure(std::string filename);

        private:
            const std::string delimiter = ",";
            const std::string moleculeDelimiter = "#";
        
    };
    
} 


#endif // !IDSIM_FILEIO_MOLECULARSTRUCTUREREADER
