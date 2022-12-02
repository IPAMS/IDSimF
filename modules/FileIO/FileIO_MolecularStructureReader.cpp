// /***************************
//  Ion Dynamics Simulation Framework (IDSimF)

//  Copyright 2020 - Physical and Theoretical Chemistry /
//  Institute of Pure and Applied Mass Spectrometry
//  of the University of Wuppertal, Germany

//  IDSimF is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.

//  IDSimF is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.

//  You should have received a copy of the GNU General Public License
//  along with IDSimF.  If not, see <https://www.gnu.org/licenses/>.
//  ****************************/

// #include "FileIO_MolecularStructureReader.hpp"
// #include "Core_constants.hpp"
// #include <fstream>
// #include <vector>
// #include <algorithm>
// #include <exception>
// #include <iostream>


// std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> FileIO::MolecularStructureReader::readMolecularStructure(std::string filename){

//     //open stream:
//     std::ifstream in;
//     in.open(filename);

//     if (in.good()){

//         //parse: read line by line, tokenize the lines
//         std::string line;
//         size_t pos = 0;
//         std::string lastMolecule;
//         std::string token;
//         std::string fullLine;
//         std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection;

//         while(std::getline(in, line)) {
//             if(line[0] == '#'){
//                 // get molecule structure key for the hash map
//                 line.erase(std::remove_if(line.begin(), line.end(), [](char chr){ return chr == '#' || chr == ',';}), line.end());
//                 line.erase(std::remove(line.begin(), line.end(), '\r' ), line.end());
//                 std::shared_ptr<CollisionModel::MolecularStructure> molstrPtr = std::make_shared<CollisionModel::MolecularStructure>();
//                 molecularStructureCollection.insert({line, molstrPtr});
//                 lastMolecule = line;
//                 molecularStructureCollection.at(lastMolecule)->setName(lastMolecule);
//                 // get diamater from next line
//                 if(std::getline(in, line)){
//                     line.erase(std::remove_if(line.begin(), line.end(), [](char chr){ return chr == ',';}), line.end());
//                     molecularStructureCollection.at(lastMolecule)->setDiameter( std::strtod(line.c_str(), nullptr) * 1E-10);
//                 }

//             }else{

//                 fullLine = line;
//                 std::string atomType;
//                 std::vector<double> atomValues = std::vector<double>();

//                 //find a token terminated by delimiter, extract it and delete it from the read line
//                 while ((pos = line.find(delimiter)) != std::string::npos) {
//                     token = line.substr(0, pos);
//                     line.erase(0, pos + delimiter.length());
//                     //line.erase(std::remove(line.begin(), line.end(), '\r' ), line.end());
//                     auto occurenceDelimiter = std::count(line.begin(), line.end(), ',');
//                     // AtomType is not a double so save this seperatly
//                     if(occurenceDelimiter == 7){
//                         atomType = token;
//                     }else{
//                         atomValues.push_back(std::strtod(token.c_str(), nullptr));
//                     }
//                 }

//                 // one token left over
//                 if (line.length() != 0) {
//                     atomValues.push_back(std::strtod(line.c_str(), nullptr));
//                 }

//                 // eight double values necessary to initialize
//                 if (atomValues.size() == 8) {
//                     std::shared_ptr<CollisionModel::Atom> atm = std::make_shared<CollisionModel::Atom>(
//                                 Core::Vector(atomValues.at(0)*1E-10, atomValues.at(1)*1E-10, atomValues.at(2)*1E-10), //location in angström -> m
//                                 atomValues.at(3), // mass in amu
//                                 atomValues.at(4), // charge in e
//                                 atomValues.at(5), // part. charge in e
//                                 CollisionModel::Atom::from_string(atomType), // Atom type
//                                 atomValues.at(6) * 1E-10,  // sigma in angström -> m
//                                 atomValues.at(7) * 1E3 / Core::N_AVOGADRO  // epsilon in kJ/mol -> J
//                     );
//                     // Insert atom to the correct molecule structure
//                     auto it = molecularStructureCollection.find(lastMolecule);
//                     if (it != molecularStructureCollection.end()) {
//                         it->second->addAtom(atm);
//                     }

//                 } else {
//                     throw  std::runtime_error("wrong number of columns in line: " + fullLine);
//                 }
//             }
//         }
//         return molecularStructureCollection;
        
//     }
//     else{
//         throw  std::runtime_error("file not found");
//     }



// }

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
#include "Core_constants.hpp"
#include <fstream>
#include <vector>
#include <algorithm>
#include <exception>
#include <iostream>
#include <json.h>


std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> FileIO::MolecularStructureReader::readMolecularStructure(std::string filename){

    //open stream:
    std::ifstream in;
    in.open(filename);

    if (in.good()){
        Json::Value confRoot;
        in >> confRoot;
        in.close();
        
        std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection;

        for (Json::Value::ArrayIndex i = 0; i != confRoot.size(); i++){
            std::shared_ptr<CollisionModel::MolecularStructure> molstrPtr = 
                                    std::make_shared<CollisionModel::MolecularStructure>();

            std::string name = confRoot[i]["name"].asString();
            double diameter = confRoot[i]["diameter"].asDouble()*1e-10;
    
            molecularStructureCollection.insert({name, molstrPtr});
            molecularStructureCollection.at(name)->setName(name);
            molecularStructureCollection.at(name)->setDiameter(diameter);
            
            for (Json::Value::ArrayIndex j = 0; j != confRoot[i]["atoms"].size(); j++){
                std::shared_ptr<CollisionModel::Atom> atm = std::make_shared<CollisionModel::Atom>(
                    Core::Vector(confRoot[i]["atoms"][j]["posx"].asDouble()*1E-10, 
                        confRoot[i]["atoms"][j]["posy"].asDouble()*1E-10, 
                        confRoot[i]["atoms"][j]["posz"].asDouble()*1E-10), //location in angström -> m
                    confRoot[i]["atoms"][j]["mass"].asDouble(), // mass in amu
                    confRoot[i]["atoms"][j]["charge"].asDouble(), // charge in e
                    confRoot[i]["atoms"][j]["partCharge"].asDouble(), // part. charge in e
                    CollisionModel::Atom::from_string(confRoot[i]["atoms"][j]["type"].asString()), // Atom type
                    confRoot[i]["atoms"][j]["LJsigma"].asDouble() * 1E-10,  // sigma in angström -> m
                    confRoot[i]["atoms"][j]["LJeps"].asDouble() * 1E3 / Core::N_AVOGADRO  // epsilon in kJ/mol -> J
                );

                auto it = molecularStructureCollection.find(name);
                if (it != molecularStructureCollection.end()) {
                    it->second->addAtom(atm);
                }
            }
        }
        
        return molecularStructureCollection;
        
    }
    else{
        throw  std::runtime_error("file not found");
    }

}