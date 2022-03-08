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
 test_exafmmtSolver.cpp

 Basic testing of ExaFMM-t based FMM solver

 ****************************/
#include "ExaFMMt_fmmSolver.hpp"
#include "Core_vector.hpp"
#include "Core_particle.hpp"
#include "PSim_boxStartZone.hpp"
#include "SC_fullSumSolver.hpp"
#include "catch.hpp"
#include "test_util.hpp"
#include "test_particleInit.hpp"
#include <iostream>

TEST_CASE( "Test particle/particle interaction calculation with ExaFMMt", "[ExaFMMt]"){

    ExaFMMt::FMMSolver exafmmSolver;
    SpaceCharge::FullSumSolver fullSumSolver;

    SECTION("Force in random box should be approximately the same as with full sum"){
        std::size_t nIons = 5000;
        Core::Vector boxSize(0.002, 0.002, 0.002);
        ParticleSimulation::BoxStartZone startZone(boxSize);
        std::vector<std::unique_ptr<Core::Particle>> ions= startZone.getRandomParticlesInStartZone(nIons, 1);

        for (std::size_t i=0; i<nIons; i++){
            exafmmSolver.insertParticle((*ions[i]), i+1);
            fullSumSolver.insertParticle((*ions[i]), i+1);
        }

        CHECK(exafmmSolver.getNumberOfParticles() == nIons);
        CHECK(fullSumSolver.getNumberOfParticles() == nIons);



        exafmmSolver.computeChargeDistribution();

        Core::Vector force1 = exafmmSolver.getEFieldFromSpaceCharge(*ions[1]);
        Core::Vector fullSumForce1 = fullSumSolver.getEFieldFromSpaceCharge(*ions[1]);
        CHECK( ((force1-fullSumForce1).magnitude() / force1.magnitude()) < 1e-2);
        std::cout << std::endl;
        std::cout << "x: "<<force1.x() / fullSumForce1.x() << std::endl;
        std::cout << "y: "<<force1.y() / fullSumForce1.y() << std::endl;
        std::cout << "z: "<<force1.z() / fullSumForce1.z() << std::endl;
        std::cout << "--------------------------------------------" << std::endl;

        Core::Vector force10 = exafmmSolver.getEFieldFromSpaceCharge(*ions[10]);
        Core::Vector fullSumForce10 = fullSumSolver.getEFieldFromSpaceCharge(*ions[10]);
        CHECK( ((force10-fullSumForce10).magnitude() / force10.magnitude()) < 1e-2);
        std::cout << "x: "<< force10.x() / fullSumForce10.x() << std::endl;
        std::cout << "y: "<<force10.y() / fullSumForce10.y() << std::endl;
        std::cout << "z: "<<force10.z() / fullSumForce10.z() << std::endl;
        std::cout << "--------------------------------------------" << std::endl;

        //CHECK( ((forceExa0-forceFullSum0).magnitude() / forceExa0.magnitude()) < 1e-2);
    }

    SECTION( "Test force calculation with a large number of particles in a latticed cube"){
        std::size_t nPerDirection = 15;
        auto ions = getIonsInLattice(nPerDirection);

        SpaceCharge::FullSumSolver fullSumSolver;
        std::size_t i = 0;
        for (auto& ion: ions){
            fullSumSolver.insertParticle(*ion, i);
            exafmmSolver.insertParticle(*ion, i);
            ++i;
        }

        exafmmSolver.computeChargeDistribution();
        CHECK(exafmmSolver.getNumberOfParticles() == nPerDirection * nPerDirection * nPerDirection);

        Core::Vector force1 = exafmmSolver.getEFieldFromSpaceCharge(*ions[1]);
        Core::Vector fullSumForce1 = fullSumSolver.getEFieldFromSpaceCharge(*ions[1]);
        CHECK( ((force1-fullSumForce1).magnitude() / force1.magnitude()) < 1e-2);
        std::cout << std::endl;
        std::cout << "x: "<<force1.x() / fullSumForce1.x() << std::endl;
        std::cout << "y: "<<force1.y() / fullSumForce1.y() << std::endl;
        std::cout << "z: "<<force1.z() / fullSumForce1.z() << std::endl;
        std::cout << "--------------------------------------------" << std::endl;

        Core::Vector force10 = exafmmSolver.getEFieldFromSpaceCharge(*ions[10]);
        Core::Vector fullSumForce10 = fullSumSolver.getEFieldFromSpaceCharge(*ions[10]);
        CHECK( ((force10-fullSumForce10).magnitude() / force10.magnitude()) < 1e-2);
        std::cout << "x: "<< force10.x() / fullSumForce10.x() << std::endl;
        std::cout << "y: "<<force10.y() / fullSumForce10.y() << std::endl;
        std::cout << "z: "<<force10.z() / fullSumForce10.z() << std::endl;
        std::cout << "--------------------------------------------" << std::endl;

        Core::Vector force100 = exafmmSolver.getEFieldFromSpaceCharge(*ions[100]);
        Core::Vector fullSumForce100 = fullSumSolver.getEFieldFromSpaceCharge(*ions[100]);
        CHECK( ((force100-fullSumForce100).magnitude() / force100.magnitude()) < 1e-2);
        std::cout << "x: "<< force100.x() / fullSumForce100.x() << std::endl;
        std::cout << "y: "<<force100.y() / fullSumForce100.y() << std::endl;
        std::cout << "z: "<<force100.z() / fullSumForce100.z() << std::endl;
        std::cout << "--------------------------------------------" << std::endl;

        /*CHECK(vectorApproxCompare(
                force1,
                fullSumForce1) == vectorsApproxEqual);*/

        //search for ion:
        double x_correct = fullSumForce1.x();

        for (auto& ion: ions){
            Core::Vector sField = exafmmSolver.getEFieldFromSpaceCharge(*ion);

            if ( std::abs(sField.x() - x_correct) < 3e-9 ){
                std::cout << fullSumForce1 << std::endl;
                std::cout << sField << std::endl;
                std::cout << "---------" << std::endl;
            }
        }

        //Core::Vector force10 = exafmmSolver.getEFieldFromSpaceCharge(*ions[10]);
        //Core::Vector fullSumForce10 = fullSumSolver.getEFieldFromSpaceCharge(*ions[10]);
        //CHECK( ((force10-fullSumForce10).magnitude() / force10.magnitude()) < 1.5e-2);
    }

}