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
 test_reactiveParticle.cpp

 Testing of reactive particle in RS

 ****************************/

#include "RS_Substance.hpp"
#include "RS_ReactiveParticle.hpp"
#include "Core_constants.hpp"
#include "catch.hpp"
#include "test_util.hpp"

TEST_CASE("Test basic instantiation and setting of species in reactive particle", "[RS][Reactive Particle]") {

    RS::Substance sub1("Testsubstance",RS::Substance::substanceType::discrete);
    RS::ReactiveParticle rp = RS::ReactiveParticle(&sub1);

    CHECK(isExactDoubleEqual(rp.getCharge(), 0.0));

    RS::Substance sub2("Testsubstance2",RS::Substance::substanceType::discrete);
    sub2.collisionDiameter(100.0);
    sub2.mass(10);
    sub2.mobility(0.1);
    sub2.charge(20);

    rp.setSpecies(&sub2);

    CHECK(isExactDoubleEqual(rp.getDiameter(), 100.0));
    CHECK(isExactDoubleEqual(rp.getMass(), 10.0*Core::AMU_TO_KG));
    CHECK(isExactDoubleEqual(rp.getMobility(), 0.1));
    CHECK(isExactDoubleEqual(rp.getCharge(), 20.0*Core::ELEMENTARY_CHARGE));

    RS::ReactiveParticle rp2 = RS::ReactiveParticle(&sub2);

    CHECK(isExactDoubleEqual(rp2.getDiameter(), 100.0));
    CHECK(isExactDoubleEqual(rp2.getMass(), 10.0*Core::AMU_TO_KG));
    CHECK(isExactDoubleEqual(rp2.getMobility(), 0.1));
    CHECK(isExactDoubleEqual(rp2.getCharge(), 20.0*Core::ELEMENTARY_CHARGE));
}