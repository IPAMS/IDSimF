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

#include "PSim_ionCloudReader.hpp"
#include "BTree_particle.hpp"
#include <fstream>

/**
 * Reads particles from an ion cloud definition file
 * @param filename a filename to read from
 * @return Vector with unique pointers to the read and newly created particles
 * @throws IonCloudFileException if the ion cloud file is incorrect
 */
std::vector<std::unique_ptr<BTree::Particle>> FileIO::IonCloudReader::readIonCloud(std::string filename) {

    //open stream:
    std::ifstream in;
    in.open(filename);

    if (in.good()){

        //parse: read line by line, tokenize the lines
        std::string line;
        std::string fullLine;

        size_t pos = 0;
        std::string token;

        std::vector<std::unique_ptr<BTree::Particle>> ionCloud;
        while(std::getline(in, line)) {
            //std::cout << line <<std::endl;
            if (line[0] != '#') {
                std::vector<double> pVal = std::vector<double>();
                fullLine = line;

                //find a token terminated by delimiter_, extract it and delete it from the read line
                while ((pos = line.find(delimiter_)) != std::string::npos) {
                    token = line.substr(0, pos);
                    line.erase(0, pos + delimiter_.length());
                    pVal.push_back(std::strtod(token.c_str(), nullptr));
                }

                //if the line was not terminated by delimiter_, a token is left:
                if (line.length() != 0) {
                    pVal.push_back(std::strtod(line.c_str(), nullptr));
                }

                if (pVal.size() == 10) {
                    ionCloud.push_back(std::make_unique<BTree::Particle>(
                                Core::Vector(pVal[0], pVal[1], pVal[2]), //location
                                Core::Vector(pVal[3], pVal[4], pVal[5]), //velocity
                                pVal[6], //charge (in elementary charges)
                                pVal[7], //mass (in amu)
                                pVal[8]*1e-10, //collision diameter (in angström)
                                pVal[9]  //time of birth
                            ));
                } else {
                    throw (FileIO::IonCloudFileException("wrong number of columns in line: " + fullLine));
                }
            }
        }
        return(ionCloud);
    }
    else{
        throw(FileIO::IonCloudFileException("file not found"));
    }
}
