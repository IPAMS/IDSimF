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

#include "FileIO_CSVReader.hpp"
#include <fstream>

std::vector<double> FileIO::CSVReader::readCSVFile(std::string filename) {

    std::ifstream in;
    in.open(filename);

    if (in.good()){
        std::string line;
        std::vector<double> result = std::vector<double>();

        while (std::getline(in, line)){
            if (line[0] != '#'){
                result.push_back(std::stod(line, nullptr));
            }
        }
        return result;
    }
    else {
        throw(FileIO::CSVReaderException("file not found"));
    }
}