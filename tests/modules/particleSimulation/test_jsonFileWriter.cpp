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
 test_jsonFileWriter.cpp

 Testing of JSON trajectory file writer

 ****************************/

#include "BTree_particle.hpp"
#include "BTree_tree.hpp"
#include "PSim_trajectoryExplorerJSONwriter.hpp"
#include "json.h"
#include <vector>
#include <iostream>
#include "catch.hpp"

Json::Value importParticleTrajectoryFile(std::string filename){
    std::ifstream particleFile;
    particleFile.open(filename);
    Json::Value root;
    particleFile>>root;
    particleFile.close();
    return root;
}

TEST_CASE("Test JSON trajectory file writer", "[ParticleSimulation][JSONTrajectoryWriter][file writers]") {

    SECTION("Json trajectory file writer should write a correct file without additional parameters") {
        std::vector<BTree::Particle*> particles;

        //create and add particles:
        BTree::Particle testIon1(Core::Vector(1.2,1.2,1.2),2.0);
        BTree::Particle testIon2(Core::Vector(1.6,1.2,1.2),2.0);
        BTree::Particle testIon3(Core::Vector(1.65,1.2,1.2),2.0);
        BTree::Particle testIon4(Core::Vector(1.66,1.2,1.2),2.0);
        BTree::Particle testIon5(Core::Vector(1.665,1.2,1.2),2.0);


        particles.push_back(&testIon1);
        particles.push_back(&testIon2);
        particles.push_back(&testIon3);
        particles.push_back(&testIon4);
        particles.push_back(&testIon5);

        // testWriter is instantiated explicitly via "new" and explicitly deleted to check if file is correctly
        // closed in the destructor of the file writer
        std::string filename("trajectoryExplorerTest.json");
        auto testWriter = new ParticleSimulation::TrajectoryExplorerJSONwriter(filename);
        testWriter->writeTimestep(particles, 0.0,false);
        testWriter->writeTimestep(particles, 0.1,false);
        testWriter->writeTimestep(particles, 0.2,true);
        delete(testWriter);

        //read back generated json file:
        Json::Value particleRoot = importParticleTrajectoryFile(filename);
        Json::Value steps = particleRoot.get("steps", 0);
        Json::Value step_2 = steps[1];

        CHECK(step_2["time"].asDouble() == 0.1);
        Json::Value step_2_pos = step_2["ions"];
        CHECK(step_2_pos[1][1].asDouble()  == 1.2);
    }

    SECTION("Json trajectory file writer should write a correct file without additional parameters") {

        std::string key_testKey = "keyTestKey";

        //Test file writer:
        std::vector<BTree::Particle*> particles;

        BTree::Particle testIon1(Core::Vector(1.2,1.2,1.2),2.0);
        BTree::Particle testIon2(Core::Vector(1.6,1.2,1.2),2.0);
        BTree::Particle testIon3(Core::Vector(1.65,1.2,1.2),2.0);

        particles.push_back(&testIon1);
        particles.push_back(&testIon2);
        particles.push_back(&testIon3);

        //root->printNode(0);
        testIon1.setFloatAttribute(key_testKey, 10.1);
        testIon2.setFloatAttribute(key_testKey, 20.2);
        testIon3.setFloatAttribute(key_testKey, 30.3);


        ParticleSimulation::partAttribTransformFctType additionalParameterTransformFct =
                [=](BTree::Particle* particle) -> std::vector<double> {
                    std::vector<double> result = {particle->getFloatAttribute(key_testKey)};
                    return result;
                };

        // testWriter is instantiated explicitly via "new" and explicitly deleted to check if file is correctly
        // closed in the destructor of the file writer
        ParticleSimulation::TrajectoryExplorerJSONwriter* testWriter =
                new ParticleSimulation::TrajectoryExplorerJSONwriter("trajectoryExplorerTest_additionalParam.json");
        std::vector<int> externalKeys = {1,2,3};

        int nSteps = 5;
        bool flagLast = false;
        for (int i=0; i<nSteps; i++){
            if (i== nSteps-1){
                flagLast = true;
            }

            testWriter->writeTimestep(
                    particles,
                    additionalParameterTransformFct,
                    std::vector<std::pair<std::string,double>> {
                            std::pair<std::string,double> {"aux_param_1",i*10},
                            std::pair<std::string,double> {"aux_param_2",i*25.1},
                    },
                    i*0.1,
                    flagLast);
        }
        delete(testWriter);

        CHECK(testIon1.getLocation() != testIon2.getLocation());

        Json::Value particleRoot = importParticleTrajectoryFile("trajectoryExplorerTest_additionalParam.json");
        Json::Value steps = particleRoot.get("steps", 0);
        Json::Value step_4 = steps[4];

        CHECK(step_4["time"].asDouble() == 0.4);
        CHECK(step_4["aux_param_2"].asDouble() == 100.4);

        Json::Value step_4_2 = step_4["ions"][2];

        CHECK(step_4_2[0][0].asDouble()  == 1.65);
        CHECK(step_4_2[1].asDouble()  == 30.3);
    }
}

