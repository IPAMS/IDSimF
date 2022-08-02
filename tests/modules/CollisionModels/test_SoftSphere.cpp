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
 test_SoftSphere.cpp

 Testing of soft sphere collision model

 ****************************/

#include "CollisionModel_SoftSphere.hpp"
#include "Core_constants.hpp"
#include "Core_randomGenerators.hpp"
#include "Core_particle.hpp"
#include "catch.hpp"


TEST_CASE("Basic test Soft Sphere model", "[CollisionModels][SoftSphereModel]") {
    //Set the global random generator to a test random number generator to make the test experiment fully deterministic:
    Core::globalRandomGeneratorPool = std::make_unique<Core::TestRandomGeneratorPool>();

    double diameterHe = CollisionModel::SoftSphereModel::DIAMETER_HE;
    Core::Particle ion;
    ion.setDiameter(CollisionModel::SoftSphereModel::DIAMETER_HE);
    ion.setVelocity(Core::Vector(100,0,0));
    ion.setMassAMU(28.0);
    int n = 100;

    SECTION("Test without maxwellian approximation"){
        CollisionModel::SoftSphereModel vss = CollisionModel::SoftSphereModel(
                1.0,298,4.0,diameterHe,false);

        for (int i =0; i<n; i++){
            vss.modifyVelocity(ion, 2e-7);
        }

        Core::Vector ionVelo = ion.getVelocity();

        //CHECK(Approx(ionVelo.x()) ==  -27.2881052);
        //CHECK(Approx(ionVelo.y()) ==  32.5138876);
        //CHECK(Approx(ionVelo.z()) ==  -128.48497);
    }
}

