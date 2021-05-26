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
 test_simulation.cpp

 Testing of chemical reaction simulations in RS

 ****************************/

#include "catch.hpp"
#include "RS_Simulation.hpp"
#include "RS_SimulationConfiguration.hpp"
#include "RS_ConfigFileParser.hpp"
#include "RS_AbstractReaction.hpp"
#include "RS_StaticReaction.hpp"
#include "Core_randomGenerators.hpp"
#include <memory>
#include <map>
#include <utility>
#include <vector>

using sMap = std::map<RS::Substance*,int>;
using sPair= sMap::value_type;


TEST_CASE( "Test basic RS simulation semantics", "[RS][Simulation]") {

    //create a simulation and add some particles:
    RS::ConfigFileParser confParser;
    RS::Simulation sim = RS::Simulation(confParser.getTestConfigWaterClusters());
    RS::Substance Cluster_1 = RS::Substance("[H3O]+",RS::Substance::substanceType::discrete);
    RS::Substance Cluster_2 = RS::Substance("[H2O+H3O]+",RS::Substance::substanceType::discrete);
    RS::Substance Cluster_3 = RS::Substance("[H4O2+H3O]+",RS::Substance::substanceType::discrete);

    uniqueReactivePartPtr p1 = std::make_unique<RS::ReactiveParticle>(&Cluster_1);
    uniqueReactivePartPtr p2 = std::make_unique<RS::ReactiveParticle>(&Cluster_2);
    uniqueReactivePartPtr p3 = std::make_unique<RS::ReactiveParticle>(&Cluster_3);



    SECTION( "Particles should be insertable into and retrievable from a simulation") {

        //check if parameters are set by the species at particle generation:
        Cluster_2.charge(5.0);
        Cluster_2.mobility(20);
        uniqueReactivePartPtr p2_2 = std::make_unique<RS::ReactiveParticle>(&Cluster_2);
        CHECK( Approx(p2_2->getCharge() - 5.0*Core::ELEMENTARY_CHARGE) == 0.0);

        sim.addParticle(p1.get(),1);
        sim.addParticle(p2_2.get(),2);


        //particles are initialized by reference, the internals of the simulation saves a reference to the particle:
        CHECK(p1.get() == &sim.getParticle(1));
        CHECK(p1.get() != &sim.getParticle(2));

        //FIXME: mobility is currently NOT taken from the species, (somewhat inconsistent currently)
        //Cluster_2.mobility(30);
        //p2->updateParticleParametersFromSpecies_();
        //CHECK(p2->getMobility() == 30);
        //CHECK(p1->getMobility() == 10);

        CHECK(p2_2->getSpecies() == &Cluster_2);

        p2_2->setSpecies(&Cluster_1);
        CHECK( p2_2->getSpecies() == &Cluster_1);

        //particle in the simulation was initalized by reference, thus particle 2 in the simulation is affected:
        CHECK( sim.getParticle(2).getSpecies() == &Cluster_1);
    }

    SECTION( "Particle double insertion in simulation is not possible") {

        p1->setChargeElementary(1.0);
        p2->setChargeElementary(2.0);
        p3->setChargeElementary(3.0);

        sim.addParticle(p1.get(),1);
        sim.addParticle(p2.get(),2);

        REQUIRE(p1->getSpecies() == sim.getParticle(1).getSpecies());

        bool insertSuccessful = sim.addParticle(p3.get(),1);

        REQUIRE(!insertSuccessful);
        REQUIRE(p1->getSpecies() == sim.getParticle(1).getSpecies());
        REQUIRE(p1->getCharge() == sim.getParticle(1).getCharge());
        REQUIRE(sim.discreteConcentrations().at(&Cluster_1) == 1);
        REQUIRE(sim.discreteConcentrations().count(&Cluster_3) == 0);
    }

    SECTION( "Discrete particle concentration calculation should be correct", "[RS Simulation]") {

        std::vector<uniqueReactivePartPtr> particles;

        for (int i=0; i < 100; ++i){
            auto particle = std::make_unique<RS::ReactiveParticle>(&Cluster_1);
            particle->setChargeElementary(i);
            particles.push_back(std::move(particle));
            sim.addParticle(particles[i].get(),i);
        }

        sim.addParticle(p2.get(), 200);

        REQUIRE(sim.discreteConcentrations().at(&Cluster_1) == 100);
        REQUIRE(sim.discreteConcentrations().at(&Cluster_2) == 1);

        for (int i=50; i < 60; ++i){
            sim.removeParticle(i);
        }

        REQUIRE(sim.discreteConcentrations().at(&Cluster_1) == 90);
        REQUIRE(sim.discreteConcentrations().at(&Cluster_2) == 1);
    }
}

TEST_CASE( "Test RS simulations", "[RS][Simulation]") {

    //Set the global random generator to a test random number generator to make the test experiment fully deterministic:
    Core::globalRandomGenerator = std::make_unique<Core::TestRandomGenerator>();
    RS::ConfigFileParser parser = RS::ConfigFileParser();


    SECTION( "Simple simulation should produce correct result") {
        RS::Simulation sim = RS::Simulation(parser.getTestConfigSimple());

        RS::Substance* subst_A = sim.simulationConfiguration()->getAllSubstances()[0];
        RS::Substance* subst_B = sim.simulationConfiguration()->getAllSubstances()[1];
        RS::Substance* subst_C = sim.simulationConfiguration()->getAllSubstances()[2];
        subst_B->staticConcentration(0.05);
        sim.simulationConfiguration()->updateConfiguration();

        int nParticles = 100;
        int nSteps = 20;
        std::vector<uniqueReactivePartPtr> particles;
        for (int i=0; i < nParticles; ++i) {
            uniqueReactivePartPtr particle = std::make_unique<RS::ReactiveParticle>(subst_A);
            sim.addParticle(particle.get(), i);
            particles.push_back(std::move(particle));
        }

        RS::ReactionConditions reactionConditions = RS::ReactionConditions();
        reactionConditions.temperature = 298;
        reactionConditions.electricField = 0.0;
        reactionConditions.pressure = 100000.0;

        for (int step=0; step < nSteps; ++step) {

            for (int i = 0; i < nParticles; ++i) {
                RS::ReactiveParticle particle = sim.getParticle(i);
                sim.react(i, reactionConditions, 1.0);
            }
            sim.advanceTimestep(1.0);
        }

        std::map<RS::Substance* const, int> finalConcs = sim.discreteConcentrations();
        REQUIRE(finalConcs[subst_A] == 37);
        REQUIRE(finalConcs[sim.simulationConfiguration()->getAllSubstances()[2]] == 63);
        REQUIRE(sim.getParticle(2).getSpecies() == subst_C);
        REQUIRE(sim.getParticle(1).getSpecies() == subst_A);

        //Check if simulation updates particles according to substances:
        REQUIRE(sim.getParticle(2).getMobility() == subst_C->mobility());
        REQUIRE(sim.getParticle(1).getMobility() == subst_A->mobility());

        REQUIRE(sim.getParticle(2).getDiameter() == subst_C->collisionDiameter());
        REQUIRE(sim.getParticle(1).getDiameter() == subst_A->collisionDiameter());

        REQUIRE(sim.getParticle(2).getCharge() == subst_C->charge()*Core::ELEMENTARY_CHARGE);
        REQUIRE(sim.getParticle(1).getCharge() == subst_A->charge()*Core::ELEMENTARY_CHARGE);

        REQUIRE(sim.getParticle(2).getMass() == subst_C->mass()*Core::AMU_TO_KG);
        REQUIRE(sim.getParticle(1).getMass() == subst_A->mass()*Core::AMU_TO_KG);
    }

    SECTION( "Simulation with water cluster configuration file should be correct") {

        RS::Simulation sim = RS::Simulation(parser.parseFile("RS_waterCluster_test.conf"));
        RS::SimulationConfiguration* simConf = sim.simulationConfiguration();

        RS::Substance* Cl1 = simConf->substanceByName("Cl_1");
        RS::Substance* H2O = simConf->substanceByName("H2O");
        H2O->staticConcentration(0.0);
        int nParticles = 100000;
        int nSteps = 200;
        double dt = 1.0e-4;
        std::vector<uniqueReactivePartPtr> particles;
        for (int i=0; i < nParticles; ++i) {
            uniqueReactivePartPtr particle = std::make_unique<RS::ReactiveParticle>(Cl1);
            sim.addParticle(particle.get(), i);
            particles.push_back(std::move(particle));
        }

        RS::ReactionConditions reactionConditions = RS::ReactionConditions();
        reactionConditions.temperature = 298;
        reactionConditions.electricField = 0.0;
        reactionConditions.pressure = 100000.0;

        for (int step=0; step < nSteps; ++step) {

            for (int i = 0; i < nParticles; ++i) {
                sim.react(i, reactionConditions, dt);
            }
            sim.advanceTimestep(dt);
        }

        REQUIRE(sim.timestep() == nSteps);
        REQUIRE(sim.simulationTime() == Approx(nSteps* dt));

        RS::AbstractReaction* reacCl1Forward = simConf->reaction(0);
        RS::AbstractReaction* reacCl2Forward = simConf->reaction(2);
        REQUIRE(reacCl1Forward->getLabel() == "cl1_forward");
        REQUIRE(reacCl2Forward->getLabel() == "cl2_forward");
        CHECK( ( sim.reactionEvents(reacCl1Forward) > 85000 && sim.reactionEvents(reacCl1Forward) < 86500) );
        CHECK( ( sim.reactionEvents(reacCl2Forward) > 58500 && sim.reactionEvents(reacCl2Forward) < 61000) );
    }

    SECTION("Result of collision based simulation with configuration file should be correct"){

        RS::Simulation sim(parser.parseFile("RS_collisionBasedReactions_test.conf"));
        int ts = sim.timestep();
        REQUIRE(ts == 0);

        RS::SimulationConfiguration* simConf = sim.simulationConfiguration();

        RS::Substance* Cl1 = simConf->substanceByName("Cl_1");
        RS::Substance* Cl2 = simConf->substanceByName("Cl_2");
        RS::Substance* N2 = simConf->substanceByName("N2");

        int nParticles = 100;
        std::vector<uniqueReactivePartPtr> particles;
        for (int i=0; i < nParticles; ++i) {
            uniqueReactivePartPtr particle = std::make_unique<RS::ReactiveParticle>(Cl2);
            sim.addParticle(particle.get(), i);
            particles.push_back(std::move(particle));
        }

        RS::CollisionConditions collisionConditions = {.totalCollisionEnergy = 0.0};

        bool hasReacted = sim.collisionReact(0,N2,collisionConditions);
        REQUIRE(hasReacted == false);
        REQUIRE(sim.getParticle(0).getSpecies() == Cl2);

        collisionConditions.totalCollisionEnergy = 9.9 / Core::JOULE_TO_EV;
        hasReacted = sim.collisionReact(1,N2,collisionConditions);
        REQUIRE(hasReacted == false);
        REQUIRE(sim.getParticle(1).getSpecies() == Cl2);

        RS::AbstractReaction* reacCl2Destruction = simConf->reaction(1);
        REQUIRE(reacCl2Destruction->getLabel() == "cl2_destruction");
        REQUIRE(sim.reactionEvents(reacCl2Destruction) == 0);

        collisionConditions.totalCollisionEnergy = 10.1 / Core::JOULE_TO_EV;
        hasReacted = sim.collisionReact(1,N2,collisionConditions);
        REQUIRE(hasReacted == true);
        REQUIRE(sim.getParticle(1).getSpecies() == Cl1);

        REQUIRE(sim.reactionEvents(reacCl2Destruction) == 1);
    }
}

