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

#include "RS_Simulation.hpp"
#include "RS_SimulationConfiguration.hpp"
#include "RS_ConfigFileParser.hpp"
#include "RS_AbstractReaction.hpp"
#include "Core_randomGenerators.hpp"
#include "catch.hpp"
#include "test_util.hpp"
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
        Cluster_2.lowFieldMobility(20);
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

        CHECK(p1->getSpecies() == sim.getParticle(1).getSpecies());

        bool insertSuccessful = sim.addParticle(p3.get(),1);

        CHECK(!insertSuccessful);
        CHECK(p1->getSpecies() == sim.getParticle(1).getSpecies());
        CHECK(isExactDoubleEqual(p1->getCharge(), sim.getParticle(1).getCharge()));
        CHECK(sim.discreteConcentrations().at(&Cluster_1) == 1);
        CHECK(sim.discreteConcentrations().count(&Cluster_3) == 0);
    }

    SECTION( "Discrete particle concentration calculation should be correct", "[RS Simulation]") {

        std::vector<uniqueReactivePartPtr> particles;

        for (unsigned int i=0; i < 100; ++i){
            auto particle = std::make_unique<RS::ReactiveParticle>(&Cluster_1);
            particle->setChargeElementary(i);
            particles.push_back(std::move(particle));
            sim.addParticle(particles[i].get(), i);
        }

        sim.addParticle(p2.get(), 200);

        CHECK(sim.discreteConcentrations().at(&Cluster_1) == 100);
        CHECK(sim.discreteConcentrations().at(&Cluster_2) == 1);

        for (std::size_t i=50; i < 60; ++i){
            sim.removeParticle(i);
        }

        CHECK(sim.discreteConcentrations().at(&Cluster_1) == 90);
        CHECK(sim.discreteConcentrations().at(&Cluster_2) == 1);
    }
}

TEST_CASE( "Test RS simulations", "[RS][Simulation]") {

    //Set the global random generator to a test random number generator to make the test experiment fully deterministic:
    Core::globalRandomGeneratorPool = std::make_unique<Core::TestRandomGeneratorPool>();
    RS::ConfigFileParser parser = RS::ConfigFileParser();


    SECTION( "Simple explicit simulation should produce correct result") {
        RS::Simulation sim = RS::Simulation(parser.getTestConfigSimple());

        RS::Substance* subst_A = sim.simulationConfiguration()->getAllSubstances()[0];
        RS::Substance* subst_B = sim.simulationConfiguration()->getAllSubstances()[1];
        RS::Substance* subst_C = sim.simulationConfiguration()->getAllSubstances()[2];
        subst_B->staticConcentration(0.05);
        sim.simulationConfiguration()->updateConfiguration();

        std::size_t nParticles = 100;
        int nSteps = 20;
        std::vector<uniqueReactivePartPtr> particles;
        for (std::size_t i=0; i < nParticles; ++i) {
            uniqueReactivePartPtr particle = std::make_unique<RS::ReactiveParticle>(subst_A);
            sim.addParticle(particle.get(), i);
            particles.push_back(std::move(particle));
        }

        RS::ReactionConditions reactionConditions = RS::ReactionConditions();
        reactionConditions.temperature = 298;
        reactionConditions.electricField = 0.0;
        reactionConditions.pressure = 100000.0;

        for (int step=0; step < nSteps; ++step) {
            for (std::size_t i = 0; i < nParticles; ++i) {
                RS::ReactiveParticle particle = sim.getParticle(i);
                sim.react(i, reactionConditions, 1.0);
            }
            sim.advanceTimestep(1.0);
        }

        std::map<RS::Substance* const, int> finalConcs = sim.discreteConcentrations();
        CHECK(finalConcs[subst_A] == 37);
        CHECK(finalConcs[sim.simulationConfiguration()->getAllSubstances()[2]] == 63);
        CHECK(sim.getParticle(2).getSpecies() == subst_C);
        CHECK(sim.getParticle(1).getSpecies() == subst_A);

        //Check if simulation updates particles according to substances:
        CHECK(isExactDoubleEqual(sim.getParticle(2).getMobility(), subst_C->lowFieldMobility()));
        CHECK(isExactDoubleEqual(sim.getParticle(1).getMobility(), subst_A->lowFieldMobility()));

        CHECK(isExactDoubleEqual(sim.getParticle(2).getDiameter(), subst_C->collisionDiameter()));
        CHECK(isExactDoubleEqual(sim.getParticle(1).getDiameter(), subst_A->collisionDiameter()));

        CHECK(isExactDoubleEqual(sim.getParticle(2).getCharge(), subst_C->charge()*Core::ELEMENTARY_CHARGE));
        CHECK(isExactDoubleEqual(sim.getParticle(1).getCharge(), subst_A->charge()*Core::ELEMENTARY_CHARGE));

        CHECK(isExactDoubleEqual(sim.getParticle(2).getMass(), subst_C->mass()*Core::AMU_TO_KG));
        CHECK(isExactDoubleEqual(sim.getParticle(1).getMass(), subst_A->mass()*Core::AMU_TO_KG));
    }

    SECTION( "Parallelized simulation with water clusters, static reaction conditions and post reaction function should be correct") {

        //switch back to real random generator, since the test generator induced artifacts with multithreading
        Core::globalRandomGeneratorPool = std::make_unique<Core::RandomGeneratorPool>();

        RS::Simulation sim = RS::Simulation(parser.parseFile("RS_waterCluster_test.conf"));
        RS::SimulationConfiguration* simConf = sim.simulationConfiguration();

        RS::Substance* Cl1 = simConf->substanceByName("Cl_1");

        std::size_t nParticles = 100000;
        int nSteps = 500;
        double dt = 6.0e-4;
        std::vector<uniqueReactivePartPtr> particles;
        for (std::size_t i=0; i < nParticles; ++i) {
            uniqueReactivePartPtr particle = std::make_unique<RS::ReactiveParticle>(Cl1);
            particle->setIntegerAttribute("reacted", 0);
            sim.addParticle(particle.get(), i);
            particles.push_back(std::move(particle));
        }

        RS::ReactionConditions reactionConditions = RS::ReactionConditions();
        reactionConditions.temperature = 298;
        reactionConditions.electricField = 0.0;
        reactionConditions.pressure = 100000.0;

        // define a post reaction function:
        long totalCountedWithFunction = 0;
        auto particlesReactedFct = [&totalCountedWithFunction](RS::ReactiveParticle* particle){
            //we had an reaction event: Count it (access to total counted value has to be synchronized)
            #pragma omp atomic
            totalCountedWithFunction++;

            //access particles
            int reactedOld = particle->getIntegerAttribute("reacted");
            particle->setIntegerAttribute("reacted", reactedOld + 1);
        };

        for (int step=0; step < nSteps; ++step) {
            sim.performTimestep(reactionConditions, dt, particlesReactedFct);
            sim.advanceTimestep(dt);
        }

        CHECK(sim.timestep() == nSteps);
        CHECK(sim.simulationTime() == Approx(nSteps* dt));

        RS::AbstractReaction* reacCl3Forward = simConf->reaction(4);
        RS::AbstractReaction* reacCl4Forward = simConf->reaction(6);
        CHECK(reacCl3Forward->getLabel() == "cl3_forward");
        CHECK(reacCl4Forward->getLabel() == "cl4_forward");

        sim.printConcentrations();
        sim.printReactionStatistics();

        CHECK( ( sim.reactionEvents(reacCl3Forward) > 383000 && sim.reactionEvents(reacCl3Forward) < 386000) );
        CHECK( ( sim.reactionEvents(reacCl4Forward) > 3022000 && sim.reactionEvents(reacCl4Forward) < 3029000) );
        CHECK( sim.totalReactionEvents() > 0);
        CHECK( sim.totalReactionEvents()  == totalCountedWithFunction);
        CHECK( particles[0]->getIntegerAttribute("reacted") > 0);
        CHECK( particles[1]->getIntegerAttribute("reacted") > 0);
    }

    SECTION( "Parallelized simulation with water clusters and individual reaction conditions should be correct") {
        RS::Simulation sim = RS::Simulation(parser.parseFile("RS_waterCluster_test_temperatureDependent.conf"));
        RS::SimulationConfiguration* simConf = sim.simulationConfiguration();

        RS::Substance* Cl6 = simConf->substanceByName("Cl_6");
        RS::Substance* Cl1 = simConf->substanceByName("Cl_1");

        std::size_t nParticles = 10000;
        int nSteps = 200;
        double dt = 1.0e-4;
        std::vector<uniqueReactivePartPtr> particles;
        for (std::size_t i=0; i < nParticles; ++i) {
            uniqueReactivePartPtr particle = std::make_unique<RS::ReactiveParticle>(Cl6);
            particle->setFloatAttribute("temperature", i*0.1);
            sim.addParticle(particle.get(), i);
            particles.push_back(std::move(particle));
        }

        auto reactionConditionsFct = [](RS::ReactiveParticle* particle, double /*time*/)->RS::ReactionConditions{
            RS::ReactionConditions reactionConditions = RS::ReactionConditions();
            reactionConditions.temperature = particle->getFloatAttribute("temperature");
            reactionConditions.electricField = 0.0;
            reactionConditions.pressure = 100000.0;

            return reactionConditions;
        };

        for (int step=0; step < nSteps; ++step) {
            sim.performTimestep(reactionConditionsFct, dt);
            sim.advanceTimestep(dt);
        }

        CHECK(sim.timestep() == nSteps);
        CHECK(sim.simulationTime() == Approx(nSteps* dt));

        auto concs = sim.discreteConcentrations();
        for (auto& co: concs){
            std::cout<<" "<< co.first->name() << " " << co.second <<std::endl;
        }

        //CHECK( ( sim.reactionEvents(reacCl1Forward) > 86000 && sim.reactionEvents(reacCl1Forward) < 87000) );
        //CHECK( ( sim.reactionEvents(reacCl2Forward) > 58500 && sim.reactionEvents(reacCl2Forward) < 61000) );
        CHECK( sim.totalReactionEvents() > 0);
        CHECK( particles[0]->getSpecies() == Cl6);
        CHECK( particles[nParticles-1]->getSpecies() == Cl1);
    }

    SECTION("Result of collision based reaction events with configuration file should be correct"){

        RS::Simulation sim(parser.parseFile("RS_collisionBasedReactions_test.conf"));
        int ts = sim.timestep();
        CHECK(ts == 0);

        RS::SimulationConfiguration* simConf = sim.simulationConfiguration();

        RS::Substance* Cl1 = simConf->substanceByName("Cl_1");
        RS::Substance* Cl2 = simConf->substanceByName("Cl_2");
        RS::Substance* N2 = simConf->substanceByName("N2");

        std::size_t nParticles = 100;
        std::vector<uniqueReactivePartPtr> particles;
        for (std::size_t i=0; i < nParticles; ++i) {
            uniqueReactivePartPtr particle = std::make_unique<RS::ReactiveParticle>(Cl2);
            sim.addParticle(particle.get(), i);
            particles.push_back(std::move(particle));
        }

        RS::CollisionConditions collisionConditions;
        collisionConditions.totalCollisionEnergy = 0.0;

        bool hasReacted = sim.collisionReact(0,N2,collisionConditions);
        CHECK(hasReacted == false);
        CHECK(sim.getParticle(0).getSpecies() == Cl2);

        collisionConditions.totalCollisionEnergy = 9.9 / Core::JOULE_TO_EV;
        hasReacted = sim.collisionReact(1,N2,collisionConditions);
        CHECK(hasReacted == false);
        CHECK(sim.getParticle(1).getSpecies() == Cl2);

        RS::AbstractReaction* reacCl2Destruction = simConf->reaction(1);
        CHECK(reacCl2Destruction->getLabel() == "cl2_destruction");
        CHECK(sim.reactionEvents(reacCl2Destruction) == 0);

        collisionConditions.totalCollisionEnergy = 10.1 / Core::JOULE_TO_EV;
        hasReacted = sim.collisionReact(1,N2,collisionConditions);
        CHECK(hasReacted == true);
        CHECK(sim.getParticle(1).getSpecies() == Cl1);

        CHECK(sim.reactionEvents(reacCl2Destruction) == 1);
    }
}

