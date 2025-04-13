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
 idealizedMultipleCollisions.cpp

 Idealized, preconstructed molecular dynamics simulations of ion-neutral collisions for a defined set of start parameters.

 ****************************/

#include "CollisionModel_MDInteractionsExperimental.hpp"
#include "CollisionModel_Molecule.hpp"
#include "CollisionModel_Atom.hpp"
#include "Core_randomGenerators.hpp"
#include "Core_constants.hpp"
#include "Core_vector.hpp"
#include "FileIO_MolecularStructureReader.hpp"
#include "appUtils_simulationConfiguration.hpp"
#include "appUtils_logging.hpp"
#include "appUtils_stopwatch.hpp"
#include "appUtils_signalHandler.hpp"
#include "appUtils_commandlineParser.hpp"
#include "FileIO_CSVReader.hpp"
#include "CollisionModel_MDForceField_Buckingham.hpp"
#include "CollisionModel_MDForceField_LJ12_6.hpp"

#include <iostream>

int main(int argc, const char * argv[]) {
    Core::globalRandomGeneratorPool = std::make_unique<Core::RandomGeneratorPool>();

    // open configuration, parse configuration file =========================================
    AppUtils::CommandlineParser cmdLineParser(argc, argv, "idealizedMultipleCollisions", "MD Simulation of a series of ion-neutral collisions", false);
    std::string projectName = cmdLineParser.resultName();
    AppUtils::logger_ptr logger = cmdLineParser.logger();

    std::string confFileName = cmdLineParser.confFileName();
    AppUtils::simConf_ptr simConf = cmdLineParser.simulationConfiguration();


    double collisionGasPolarizability_m3 = simConf->doubleParameter("collision_gas_polarizability_m3");
    std::string collisionGasIdentifier = simConf->stringParameter("collision_gas_identifier");
    std::string particleIdentifier = simConf->stringParameter("particle_identifier");
    std::string potentialsFF = simConf->stringParameter("force_field");
    std::string outputFilename = simConf->stringParameter("output_file");


    
    // ======================================================================================

    //read position and velocity input 
    std::vector<Core::Vector> positions; 

    
    std::string startingConfiguration = simConf->pathRelativeToConfFile(simConf->stringParameter("starting_configuration"));
    FileIO::CSVReader startConfReader = FileIO::CSVReader();
    std::vector<std::vector<std::string>> startingConfigurationCollection = startConfReader.readCSVFile(startingConfiguration, ' ');
    for(auto &line : startingConfigurationCollection){
        positions.push_back({std::stod(line[0]), std::stod(line[1]), std::stod(line[2])});
    }
    //read molecular structure file
    
    std::string mdCollisionConfFile = simConf->pathRelativeToConfFile(simConf->stringParameter("md_configuration"));
    FileIO::MolecularStructureReader mdConfReader = FileIO::MolecularStructureReader();
    std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection = 
                                                                        mdConfReader.readMolecularStructure(mdCollisionConfFile);

    

    //auto forceFieldPtr = nullptr;
    CollisionModel::MDForceField_LJ12_6 forceField(collisionGasPolarizability_m3, potentialsFF);
    auto forceFieldPtr = std::make_unique<CollisionModel::MDForceField_LJ12_6>(forceField);

    CollisionModel::Molecule ion({0,0,0}, {0,0,0}, molecularStructureCollection.at(particleIdentifier));

    CollisionModel::Molecule bgGas({0,0,0}, {0,0,0}, molecularStructureCollection.at(collisionGasIdentifier));

    static std::ofstream forcesOut;
    if(!forcesOut.is_open()){
        forcesOut.open(outputFilename);
    }


    size_t samples = positions.size();
    forcesOut << "x, y, z, f_x, f_y, f_z, vdw_x, vdw_y, vdw_z, ii_x, ii_y, ii_z" << std::endl;
    for(size_t i = 0; i < samples; i++){
        bgGas.setComPos(positions[i]);
        std::vector<CollisionModel::Molecule*> moleculesPtr = {&ion, &bgGas};
        std::vector<Core::Vector> forceMolecules(moleculesPtr.size());
        Core::Vector forceVDW;
        Core::Vector forceII;
        forceFieldPtr->calculateForceFieldComponents(moleculesPtr, forceMolecules, forceVDW, forceII);
        forcesOut << positions[i].x() << "," << positions[i].y() << "," << positions[i].z() 
        << "," << forceMolecules[1].x() << "," << forceMolecules[1].y() << "," << forceMolecules[1].z()
        << "," << forceVDW.x() << "," << forceVDW.y() << "," << forceVDW.z()
        << "," << forceII.x() << "," << forceII.y() << "," << forceII.z() << std::endl;
    }
}

