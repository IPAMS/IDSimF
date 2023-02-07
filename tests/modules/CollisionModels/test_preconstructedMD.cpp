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
 test_MDInteractions.cpp

 Testing of molecular collision model with LJ-12-6 and dipole forces

 ****************************/

#include "CollisionModel_MDInteractionsPreconstructed.hpp"
#include "CollisionModel_Molecule.hpp"
#include "CollisionModel_Atom.hpp"
#include "Core_randomGenerators.hpp"
#include "Core_constants.hpp"
#include "catch.hpp"
#include "FileIO_MolecularStructureReader.hpp"
#include <iostream>


TEST_CASE("Basic test MD preconstructed", "[CollisionModels][MDInteractionsModel]") {

    Core::globalRandomGeneratorPool = std::make_unique<Core::TestRandomGeneratorPool>();

    double diameterN2 = CollisionModel::MDInteractionsModelPreconstructed::DIAMETER_N2;
    FileIO::MolecularStructureReader reader = FileIO::MolecularStructureReader();
    std::unordered_map<std::string,  std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection = 
                                                                    reader.readMolecularStructure("test_molecularstructure_reader.json");
    Core::Particle ion;
    ion.setMolecularStructure(molecularStructureCollection.at("Li+"));
    Core::RandomSource* rndSource = Core::globalRandomGeneratorPool->getThreadRandomSource();
    
    int samples = 3;
    double pi = 3.1415;
    std::vector<Core::Vector> ionVelocities;
    std::vector<Core::Vector> ionPositions;
    std::vector<Core::Vector> ionRotations;
    for(int i = 0; i < samples; i++){
        Core::Vector tmp1 = Core::Vector(rndSource->uniformRealRndValue()*100, 0.0, 0.0);
        ionVelocities.push_back(tmp1);
        ionPositions.push_back(Core::Vector(40e-10, 0.0, 0.0));
        // Core::Vector tmp2 = Core::Vector(0.0, 0.0, rndSource->uniformRealRndValue()*1*pi/2-pi/4);
        Core::Vector tmp2 = Core::Vector(0.0, pi/2, pi/2);
        ionRotations.push_back(tmp2);
    }
    
    for(size_t i = 0; i < ionVelocities.size(); i++){
        ion.setVelocity(ionVelocities[i]);
        CollisionModel::MDInteractionsModelPreconstructed mdSim = CollisionModel::MDInteractionsModelPreconstructed(2000000, 298, 28, 
                                                                                        diameterN2,
                                                                                        1.7E-30, 
                                                                                        "N2", 
                                                                                        1e-9, 
                                                                                        1E-17, 
                                                                                        3, 1, 
                                                                                        40e-10, 
                                                                                        molecularStructureCollection, 
                                                                                        ionPositions[i],
                                                                                        ionRotations[i]);

        
        mdSim.setTrajectoryWriter("MD_collisions_preconstructed_COMeq_trajectories.txt", 40e-10, 0);
        mdSim.updateModelTimestepParameters(1, 0);
        double dt = 2e-11;
        mdSim.modifyVelocity(ion, dt);
    }
    


}
