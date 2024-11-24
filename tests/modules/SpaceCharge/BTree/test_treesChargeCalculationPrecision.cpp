/***************************
Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2024 - Physical and Theoretical Chemistry /
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
 test_treesChargeCalculationPrecision.cpp

 Testing of precision / validity of force calculation with serial and parallel tree implementations

 ****************************/

#include "Core_vector.hpp"
#include "Core_particle.hpp"
#include "BTree_tree.hpp"
#include "BTree_parallelTree.hpp"
#include "PSim_boxStartZone.hpp"
#include "PSim_util.hpp"
#include "SC_fullSumSolver.hpp"
#include "test_particleInit.hpp"
#include "catch.hpp"
#include "test_util.hpp"
#include <iostream>

void testChargedGrids(Core::Vector loc_min, Core::Vector loc_max, double chargeTop, double chargeBottom, double theta, bool print) {

    BTree::Tree testTree(loc_min,loc_max);
    BTree::ParallelTree testTreeParallel(loc_min,loc_max);

    testTreeParallel.getRoot()->setTheta(theta);

    CHECK(testTree.getRoot()->getTheta() == Approx(theta));
    unsigned int nPerDirection = 30;
    auto cations = getIonsOnXYGrid(nPerDirection, 0.004, 0.001, chargeTop);
    auto anions = getIonsOnXYGrid(nPerDirection, 0.004, -0.001, chargeBottom);
    std::vector<std::unique_ptr<Core::Particle>> ions;
    SpaceCharge::FullSumSolver fullSumSolver;

    std::size_t i = 0;
    for (auto& ion: cations){
        if (print) {
            std::cout << ion->getLocation()<<" | ";
        }
        fullSumSolver.insertParticle(*ion, i);
        testTree.insertParticle(*ion, i);
        testTreeParallel.insertParticle(*ion, i);
        ions.emplace_back(std::move(ion));
        ++i;
    }
    std::cout << std::endl;

    for (auto& ion: anions){
        if (print) {
            std::cout << ion->getLocation()<<" | ";
        }
        fullSumSolver.insertParticle(*ion, i);
        testTree.insertParticle(*ion, i);
        testTreeParallel.insertParticle(*ion, i);
        ions.emplace_back(std::move(ion));
        ++i;
    }
    std::cout << std::endl;

    Core::Particle testParticle1({0.0001, 0.0005, -0.0006}, +1);
    Core::Particle testParticle2({0.0001, 0.0005, -0.0006}, -1);

    ions.emplace_back(std::make_unique<Core::Particle>(testParticle1));
    ions.emplace_back(std::make_unique<Core::Particle>(testParticle2));

    testTree.computeChargeDistribution();
    testTreeParallel.init();
    CHECK(testTree.getNumberOfParticles() == nPerDirection * nPerDirection * 2);

    std::vector<std::size_t> ionsToTest;// = {1,4,8,9};
    ionsToTest.push_back(2);
    //ionsToTest.push_back(500);
    //ionsToTest.push_back(1500);
    //ionsToTest.push_back(1300);
    //ionsToTest.push_back(1100);

    for (auto& ionToTest: ionsToTest) {
        double charge = ions[ionToTest]->getCharge();
        if (print) {
            std::cout <<"i:   "<<ionToTest<<"  location: "<<ions[ionToTest]->getLocation()<<"  charge: "<<charge<<std::endl;
        }
        Core::Vector force = testTree.getEFieldFromSpaceCharge(*ions[ionToTest])*charge;
        Core::Vector forceParallel = testTreeParallel.getEFieldFromSpaceCharge(*ions[ionToTest])*charge;
        Core::Vector fullSumForce = fullSumSolver.getEFieldFromSpaceCharge(*ions[ionToTest])*charge;

        CHECK(force.x() == Approx(forceParallel.x()));
        CHECK(force.y() == Approx(forceParallel.y()));
        CHECK(force.z() == Approx(forceParallel.z()));

        //FIXME: Add check of deviation from full force
        if (print) {
            std::cout <<"bt force:   "<<force<<"\n"<<"bt force p: "<<forceParallel<<"\n"<<"full force: "<<fullSumForce<<"\n"<<"---------------------"<<std::endl;
        }
        //CHECK( ((force-fullSumForce).magnitude() / force.magnitude()) < 1e-2);
    }
}


TEST_CASE( "Test serial tree charge distribution calculation bipolar","[Tree]") {
    Core::Vector loc_min = Core::Vector(-1000,-1000,-1000);
    Core::Vector loc_max = Core::Vector( 1000, 1000, 1000);

    SECTION( "Test force calculation with high theta"){
        testChargedGrids(loc_min, loc_max,  1, 1, 0.9, true);
        testChargedGrids(loc_min, loc_max, -1,-1, 0.9, true);
        testChargedGrids(loc_min, loc_max, -1, 1, 0.9, true);
        testChargedGrids(loc_min, loc_max,  1,-1, 0.9, true);
    }

    SECTION( "Test force calculation with positive and negative charges medium theta") {
        testChargedGrids(loc_min, loc_max,  1, 1, 0.5, true);
        testChargedGrids(loc_min, loc_max, -1,-1, 0.5, true);
        testChargedGrids(loc_min, loc_max, -1, 1, 0.5, true);
        testChargedGrids(loc_min, loc_max,  1,-1, 0.5, true);
    }

    SECTION( "Test force calculation with positive and negative charges low theta") {
        testChargedGrids(loc_min, loc_max,  1, 1, 0.1, true);
        testChargedGrids(loc_min, loc_max, -1,-1, 0.1, true);
        testChargedGrids(loc_min, loc_max, -1, 1, 0.1, true);
        testChargedGrids(loc_min, loc_max,  1,-1, 0.1, true);;
    }
}