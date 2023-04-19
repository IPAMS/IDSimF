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

            if(name == "N2"){
                Core::Vector bondVector = molstrPtr->getAtoms().at(0)->getRelativePosition() - molstrPtr->getAtoms().at(1)->getRelativePosition();
                double bondLength = sqrt(bondVector.x()*bondVector.x() + bondVector.y()*bondVector.y() + bondVector.z()*bondVector.z());

                // this assumes the first two atoms are the nitrogen atoms!! 
                molstrPtr->getAtoms().at(0)->getRelativePosition().x(bondLength/2);
                molstrPtr->getAtoms().at(0)->getRelativePosition().y(0);
                molstrPtr->getAtoms().at(0)->getRelativePosition().z(0);
                molstrPtr->getAtoms().at(1)->getRelativePosition().x(-bondLength/2);
                molstrPtr->getAtoms().at(1)->getRelativePosition().y(0);
                molstrPtr->getAtoms().at(1)->getRelativePosition().z(0);
            }
        }
        
        return molecularStructureCollection;
        
    }
    else{
        throw  std::runtime_error("file not found");
    }

}