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
 FileIO_CSVReader.hpp

 Simple reader class for csv files

 ****************************/

#ifndef IDSIM_FILEIO_CSVREADER
#define IDSIM_FILEIO_CSVREADER

#include <iostream>
#include <vector>

namespace FileIO{
class CSVReaderException : public std::runtime_error {
    public:
        explicit CSVReaderException (const std::string msg): std::runtime_error(msg){}
    };

class CSVReader {
    public:
        [[nodiscard]] std::vector<std::vector<std::string>> readCSVFile(std::string filename, char delimiter);
        [[nodiscard]] std::vector<double> extractDouble(std::vector<std::vector<std::string>> &string_vector, unsigned int column);
        [[nodiscard]] std::vector<std::string> extractString(std::vector<std::vector<std::string>> &string_vector, unsigned int column);
    };
}

#endif /* IDSIM_FILEIO_CSVREADER */