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

#include "FileIO_MolecularStructureReader.hpp"
#include "CollisionModel_MolecularStructure.hpp"
#include "Core_constants.hpp"
#include <fstream>
#include <vector>
#include <algorithm>
#include <exception>
#include <memory>


void FileIO::MolecularStructureReader::readMolecularStructure(std::string filename){

    //open stream:
    std::ifstream in;
    in.open(filename);

    if (in.good()){

        //parse: read line by line, tokenize the lines
        std::string line;
        size_t pos = 0;
        std::string lastMolecule;
        std::string token;
        std::string fullLine;

        while(std::getline(in, line)) {
            if(line[0] == '#'){
                // get molecule structure key for the hash map
                line.erase(std::remove_if(line.begin(), line.end(), [](char chr){ return chr == '#' || chr == ',';}), line.end());
                std::shared_ptr<CollisionModel::MolecularStructure> molstrPtr = std::make_shared<CollisionModel::MolecularStructure>();
                CollisionModel::MolecularStructure::molecularStructureCollection.insert({line, molstrPtr});
                lastMolecule = line;
                // get diamater from next line
                if(std::getline(in, line)){
                    line.erase(std::remove_if(line.begin(), line.end(), [](char chr){ return chr == ',';}), line.end());
                    CollisionModel::MolecularStructure::molecularStructureCollection.at(lastMolecule)->setDiameter( std::strtod(line.c_str(), nullptr) * 1E-10);
                }

            }else{

                fullLine = line;
                std::string atomType;
                std::vector<double> atomValues = std::vector<double>();

                //find a token terminated by delimiter, extract it and delete it from the read line
                while ((pos = line.find(delimiter)) != std::string::npos) {
                    token = line.substr(0, pos);
                    line.erase(0, pos + delimiter.length());
                    auto occurenceDelimiter = std::count(line.begin(), line.end(), ',');
                    // AtomType is not a double so save this seperatly
                    if(occurenceDelimiter == 7){
                        atomType = token;
                    }else{
                        atomValues.push_back(std::strtod(token.c_str(), nullptr));
                    }
                }

                // one token left over
                if (line.length() != 0) {
                    atomValues.push_back(std::strtod(line.c_str(), nullptr));
                }

                // eight double values necessary to initialize
                if (atomValues.size() == 8) {
                    std::shared_ptr<CollisionModel::Atom> atm = std::make_shared<CollisionModel::Atom>(
                                Core::Vector(atomValues.at(0)*1E-10, atomValues.at(1)*1E-10, atomValues.at(2)*1E-10), //location in angström
                                atomValues.at(3), // mass in amu
                                atomValues.at(4), // charge in e
                                atomValues.at(5), // part. charge in e
                                CollisionModel::Atom::from_string(atomType), // Atom type
                                atomValues.at(6) * 1E-10,  // sigma in angström
                                atomValues.at(7) * 1E3 / Core::N_AVOGADRO  // epsilon in kJ/mol
                    );
                    // Insert atom to the correct molecule structure
                    auto it = CollisionModel::MolecularStructure::molecularStructureCollection.find(lastMolecule);
                    if (it != CollisionModel::MolecularStructure::molecularStructureCollection.end()) {
                        it->second->addAtom(atm);
                    }

                } else {
                    throw  std::runtime_error("wrong number of columns in line: " + fullLine);
                }
            }
        }
        
    }
    else{
        throw  std::runtime_error("file not found");
    }



}