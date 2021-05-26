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
 test_reaction.cpp

 Testing of chemical reactions in RS

 ****************************/

#include "catch.hpp"
#include "RS_Substance.hpp"
#include "RS_StaticReaction.hpp"
#include "RS_StaticThermalizingReaction.hpp"
#include "RS_VantHoffReaction.hpp"
#include "RS_FieldDependentVantHoffReaction.hpp"
#include "RS_SimpleCollisionStepReaction.hpp"
#include "Core_constants.hpp"
#include "Core_randomGenerators.hpp"

#include <array>
#include <numeric>


using sMap = std::map<RS::Substance*,int>;
using sMapPtr = const sMap*;
using sPair= sMap::value_type;

TEST_CASE("Test basic reaction semantics", "[RS][Reaction]") {
    SECTION("RS Reactions should be usable as map keys") {

        std::map<RS::Substance, int> st;
        RS::Substance sub_1 = RS::Substance("AAA_testsubstance", RS::Substance::substanceType::discrete);
        RS::Substance sub_2 = RS::Substance("BBB_testsubstance", RS::Substance::substanceType::discrete);

        st[sub_1] = 1;
        CHECK(st[sub_1] == 1);
        st[sub_1] = 2;
        CHECK(st[sub_1] == 2);
        st[sub_2] = 4;
        CHECK(st[sub_2] == 4);
    }

    SECTION("RS ubstances should be usable in substance table of a RS reaction") {

        RS::Substance ed_1 = RS::Substance("educt_1",RS::Substance::substanceType::discrete);
        RS::Substance ed_2 = RS::Substance("educt_2",RS::Substance::substanceType::discrete);
        RS::Substance ed_3 = RS::Substance("educt_3",RS::Substance::substanceType::isotropic);
        RS::Substance ed_4 = RS::Substance("educt_4",RS::Substance::substanceType::isotropic);
        RS::Substance pro_1 = RS::Substance("product_1",RS::Substance::substanceType::discrete);
        RS::Substance pro_2 = RS::Substance("product_2",RS::Substance::substanceType::discrete);
        RS::Substance pro_3 = RS::Substance("product_3",RS::Substance::substanceType::discrete);
        RS::Substance pro_4 = RS::Substance("product_4",RS::Substance::substanceType::isotropic);

        sMap educts;
        sMap products;

        educts.insert(sPair(&ed_1,2));
        educts.insert(sPair(&ed_2,3));
        educts.insert(sPair(&ed_3,5));
        products.insert(sPair(&pro_1,4));
        products.insert(sPair(&pro_2,2));

        RS::StaticReaction reac = RS::StaticReaction(educts, products, 1.0, "a test reaction");

        CHECK(reac.getLabel() == "a test reaction");

        sMapPtr ret_dist_educts = reac.discreteEducts();
        sMapPtr ret_dist_products = reac.discreteProducts();

        CHECK(ret_dist_educts->at(&ed_1) == 2);
        CHECK(ret_dist_educts->count(&ed_3) == 0);
        CHECK(ret_dist_educts->count(&ed_4) == 0);

        CHECK(ret_dist_products->at(&pro_1) == 4);
        CHECK(ret_dist_products->at(&pro_2) == 2);
        CHECK(ret_dist_products->count(&pro_3) == 0);
        CHECK(ret_dist_products->count(&pro_4) == 0);
    }
}

TEST_CASE("Test chemical semantics of RS reaction types", "[RS][Reaction]") {

    // setup an educt and a product substance:
    RS::Substance ed_1 = RS::Substance("educt_1", RS::Substance::substanceType::discrete);
    RS::Substance pro_1 = RS::Substance("product_1", RS::Substance::substanceType::discrete);

    sMap educts;
    sMap products;

    educts.insert(sPair(&ed_1,0));
    products.insert(sPair(&pro_1,0));

    RS::ReactionConditions reactionConditions;
    RS::ReactiveParticle dummyParticle(&ed_1);

    SECTION("Reaction probability of static reaction should be correct") {
        RS::StaticReaction reac = RS::StaticReaction(educts, products, 1.5e10, "a test reaction");
        CHECK(reac.attemptReaction(reactionConditions, &dummyParticle, 1.0).reactionProbability == 1.5e10);
        CHECK(reac.attemptReaction(reactionConditions, &dummyParticle, 0.1).reactionProbability == 1.5e9);
    }

    SECTION("Static thermalizing reaction should calculate correct reaction probabilities and thermalize reacted particle") {

        //This test is a statistical test: use a real random number generator:
        Core::globalRandomGenerator = std::make_unique<Core::RandomGenerator>();

        pro_1.mass(200);
        RS::StaticThermalizingReaction reac = RS::StaticThermalizingReaction(educts, products, 10.0, "a test reaction");
        reactionConditions.temperature = 298;

        // test if linearized reaction probability is correct:
        RS::ReactiveParticle productParticle(&pro_1);
        CHECK(reac.attemptReaction(reactionConditions, &productParticle, 1.0).reactionProbability == 10.0);

        // test if the velocity is reinitialized thermally:
        // Generate 100000 test samples and determine mean velocity and mean velocity magnitude after collision
        int nSamples = 100000;
        std::vector<Core::Vector> velocities(nSamples);
        std::vector<double> magnitudes;

        for (auto& vel: velocities){
            productParticle.setVelocity({1000.0, 0.0, 0.0});
            reac.attemptReaction(reactionConditions, &productParticle, 1.0);
            vel = productParticle.getVelocity();
        }

        Core::Vector meanVelocity =
                std::accumulate(
                        velocities.begin(), velocities.end(), Core::Vector(0.0,0.0,0.0)) /
                        velocities.size();

        REQUIRE(meanVelocity.magnitude() < 1.0);

        std::transform(velocities.begin(), velocities.end(), std::back_inserter(magnitudes),
                [](Core::Vector vec) -> double { return vec.magnitude(); });

        double meanVelocityMagnitude =
                std::accumulate(
                        magnitudes.begin(), magnitudes.end(), 0.0) /
                        magnitudes.size();


        //calculate analytical result:
        double analyticalMeanVelocityMagnitude = std::sqrt(
                (8.0 * RS::kBoltzmann * reactionConditions.temperature) / (M_PI * productParticle.getMass() ));

        REQUIRE(meanVelocityMagnitude == Approx(analyticalMeanVelocityMagnitude).margin(2));
    }

    SECTION( "Reaction probability of a van't Hoff dependent reaction should be correct") {

        RS::VantHoffReaction reac = RS::VantHoffReaction(
                educts,products,
                132.0*1000.0,6.8408e17,2e-9,
                "a test reaction");

        reactionConditions.temperature = 298;
        double a = reac.attemptReaction(reactionConditions, &dummyParticle, 1.0).reactionProbability;
        REQUIRE(a - 2.84631e-27 < 1e-30);

        reactionConditions.temperature = 398;
        double b = reac.attemptReaction(reactionConditions, &dummyParticle, 1.0).reactionProbability;
        REQUIRE(b - 1.85175e-21 < 1e-25);
    }

    SECTION("Reaction probability of field dependent van't hoff reaction rate should be correct") {

        RS::FieldDependentVantHoffReaction reac = RS::FieldDependentVantHoffReaction(
                educts,products,
                49.0*1000.0,
                1.0189e+03,
                2.0e-9,
                2.085e-4,
                28,
                "a test reaction (water cluster 1)");
        REQUIRE(!reac.isCollisionReaction());

        reactionConditions.temperature = 298;
        reactionConditions.electricField = 0;
        reactionConditions.pressure = 101325;

        double k = reac.attemptReaction(reactionConditions, &dummyParticle, 1.0).reactionProbability;

        REQUIRE(k - 1.94347e-12 < 1e-16);

        reactionConditions.electricField = 10000;
        k = reac.attemptReaction(reactionConditions, &dummyParticle, 1.0).reactionProbability;
        REQUIRE(k - 1.9441e-12 < 1e-16);

        reactionConditions.pressure = 10132.5;
        k = reac.attemptReaction(reactionConditions, &dummyParticle, 1.0).reactionProbability;
        REQUIRE(k - 2.00727e-12 < 1e-16);
    }

    SECTION("Reaction probability of collision based step reaction should be correct") {

        RS::SimpleCollisionStepReaction reac(
                educts,products,
                1.5e-1,
                "a simple collision test reaction");

        REQUIRE(reac.isCollisionReaction());

        //reaction below activation energy should return zero:
        double totalReactionEnergy = 1.0e-1 / Core::JOULE_TO_EV;
        RS::CollisionConditions collisionConditions = {.totalCollisionEnergy = totalReactionEnergy};
        double probability = reac.attemptReaction(collisionConditions,&dummyParticle).reactionProbability;
        REQUIRE(probability ==Approx(0.0));

        //reaction above activation energy should return preexponential factor * isotropic educts concentration * dt
        totalReactionEnergy = 2.0e-1 / Core::JOULE_TO_EV;
        collisionConditions.totalCollisionEnergy = totalReactionEnergy;
        probability = reac.attemptReaction(collisionConditions, &dummyParticle).reactionProbability;
        REQUIRE(probability ==Approx(1.0));
    }
}