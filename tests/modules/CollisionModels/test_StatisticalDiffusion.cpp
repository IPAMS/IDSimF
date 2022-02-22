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
 test_StatisticalDiffusion.cpp

 Testing of statistical diffusion simulation (SDS) collsion model

 ****************************/
#include "Core_randomGenerators.hpp"
#include "BTree_particle.hpp"
#include "BTree_tree.hpp"
#include "CollisionModel_StatisticalDiffusion.hpp"
#include "CollisionStatistic_default.hpp"
#include "PSim_constants.hpp"
#include "PSim_verletIntegrator.hpp"
#include "catch.hpp"


TEST_CASE( "Test SDS collison model", "[CollisionModels][SDS]") {

    //Set the global random generator to a test random number generator to make the test experiment fully deterministic:
    Core::globalRandomGeneratorPool = std::make_unique<Core::TestRandomGeneratorPool>();

    SECTION( "SDS random jumps should be correct with default statistics") {
        CollisionModel::StatisticalDiffusionModel sds(100000,298,28,0.366 * 1.0e-9);
        BTree::Particle ion;
        ion.setMassAMU(100);
        sds.setSTPParameters(ion);
        sds.updateModelParameters(ion);
        for (int i=0; i< 100; i++){
            //Core::Vector position = ion.getLocation();
            sds.modifyPosition(ion.getLocation(),ion,1e-2);
            //ion.setLocation(position);
        }
        Core::Vector ionLoc = ion.getLocation();
        CHECK(Approx(ionLoc.x()) == 0.0451560);
        CHECK(Approx(ionLoc.y()) == 0.0164318);
        CHECK(Approx(ionLoc.z()) == 0.0204592);
    }

    SECTION( "SDS random jumps with modified statistics should be correct") {
        CollisionModel::CollisionStatistics cs("cs_icdf_2020_03_15_001_parameterTest.dat");
        CollisionModel::StatisticalDiffusionModel sds(100000,298,2,0.366 * 1.0e-9, cs);
        BTree::Particle ion;
        ion.setMassAMU(12);
        sds.setSTPParameters(ion);
        sds.updateModelParameters(ion);
        for (int i=0; i< 100; i++){
            //Core::Vector position = ion.getLocation();
            sds.modifyPosition(ion.getLocation(),ion,1e-2);
            //ion.setLocation(position);
        }
        Core::Vector ionLoc = ion.getLocation();
        CHECK(Approx(ionLoc.x()) == 0.0884540);
        CHECK(Approx(ionLoc.y()) == 0.0316142);
        REQUIRE(Approx(ionLoc.z()) == 0.0400663);
    }

    SECTION( "Calculation of mean free path and thermal velocity with varied MFP "
             "and mobility with default statistics should be correct") {

        double massGas_amu = 28.94515; //SDS default value for air
        double diameterGas = 0.366 * 1.0e-9;

        CollisionModel::StatisticalDiffusionModel sds(100000,298.15,massGas_amu,diameterGas);
        BTree::Particle ion;
        //n=6, m=91, d=0.541621(est), Ko=1.89928(est), MFPo=3.2037e-005, Vo=0.252097 [*]
        //double massIon_amu = 91;
        //double diameterIon = 0.541621 * 1.0e-9;
        double massIon_amu = massGas_amu;
        double diameterIon = diameterGas;

        ion.setDiameter(diameterIon);
        ion.setMassAMU(massIon_amu);

        sds.setSTPParameters(ion);
        REQUIRE(Approx(ion.getMeanThermalVelocitySTP()) == 447);
        REQUIRE(std::abs(ion.getMeanFreePathSTP() - 6.25376e-08) < 1e-12);
    }

    SECTION( "Simulation of stokes style damping with random jumps and default statistics should be correct") {

        //Test with verlet integration:
        double dt = 1e-4;
        unsigned int timeSteps = 50;

        double massGas_amu = 28.94515; //SDS default value for air
        double diameterGas = 0.366 * 1.0e-9;

        BTree::Particle ion;

        double massIon_amu = massGas_amu;
        double diameterIon = diameterGas;

        ion.setDiameter(diameterIon);
        ion.setMassAMU(massIon_amu);
        ion.setMobility(3.5e-4);

        std::vector<BTree::Particle*>particles= std::vector<BTree::Particle*>();
        particles.push_back(&ion);
        CollisionModel::StatisticalDiffusionModel sds(100000,298.15,massGas_amu,diameterGas);
        sds.setSTPParameters(ion);
        //sds.updateModelParameters(ion);

        double ionAcceleration = 9.6485335276142e9; //((1000V / 100mm) * elementary charge) / 100 amu = 9.64e9 m/s^2

        auto accelerationFct = [ionAcceleration] (BTree::Particle* /*particle*/, int /*particleIndex*/,
                BTree::Tree& /*tree*/, double /*time*/, int /*timestep*/){
            Core::Vector result(ionAcceleration,0,0);
            return(result);
        };

        Integration::VerletIntegrator verletIntegrator(
                particles,
                accelerationFct,
                ParticleSimulation::noFunction, ParticleSimulation::noFunction, ParticleSimulation::noFunction,
                &sds);

        verletIntegrator.run(timeSteps, dt);

        Core::Vector ionPos = ion.getLocation();

        REQUIRE(Approx(ionPos.x()).epsilon(1e-3) == 0.006103848);
        REQUIRE(Approx(ionPos.y()).epsilon(1e-3) == 0.000004531);
        REQUIRE(Approx(ionPos.z()).epsilon(1e-3) == 0.000304865);
    }
}