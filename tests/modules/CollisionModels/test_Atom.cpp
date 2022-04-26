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
 test_Atom.cpp

 Testing of Atom class for use in MD Interactions Collision Model.

 ****************************/

#include "CollisionModel_Atom.hpp"
#include "Core_constants.hpp"
#include "catch.hpp"
#include "test_util.hpp"

TEST_CASE("Basic Atom construction tests", "[CollisionModels][Atom]") {

    SECTION("Default constructor"){
        CollisionModel::Atom atm = CollisionModel::Atom();

        CHECK(atm.getRelativePosition() == Core::Vector(0.0, 0.0, 0.0));
        CHECK(isExactDoubleEqual(atm.getMass(), 0.0*Core::AMU_TO_KG));
        CHECK(isExactDoubleEqual(atm.getCharge(), 0.0*Core::ELEMENTARY_CHARGE));
        CHECK(isExactDoubleEqual(atm.getPartCharge(), 0.0*Core::ELEMENTARY_CHARGE));
        CHECK(atm.getType() == CollisionModel::Atom::AtomType::H);
        CHECK(isExactDoubleEqual(atm.getSigma(), 0.0));
        CHECK(isExactDoubleEqual(atm.getEpsilon(), 0.0));

    }

    SECTION("Constructor with position, mass and charge"){
        CollisionModel::Atom atm = CollisionModel::Atom(Core::Vector(0.5, 0.1, -0.2), 4, -1.0);

        CHECK(atm.getRelativePosition() == Core::Vector(0.5, 0.1, -0.2));
        CHECK(isExactDoubleEqual(atm.getCharge(), -1.0*Core::ELEMENTARY_CHARGE));
        CHECK(isExactDoubleEqual(atm.getMass(), 4*Core::AMU_TO_KG));
        CHECK(isExactDoubleEqual(atm.getPartCharge(), 0.0*Core::ELEMENTARY_CHARGE));
        CHECK(atm.getType() == CollisionModel::Atom::AtomType::H);
        CHECK(isExactDoubleEqual(atm.getSigma(), 0.0));
        CHECK(isExactDoubleEqual(atm.getEpsilon(), 0.0));
    }

    SECTION("Constructor with position, mass and (partial) charge"){
        CollisionModel::Atom atm = CollisionModel::Atom(Core::Vector(0.5, 0.1, -0.2), 4, -1.0, 2.0);

        CHECK(atm.getRelativePosition() == Core::Vector(0.5, 0.1, -0.2));
        CHECK(isExactDoubleEqual(atm.getCharge(), -1.0*Core::ELEMENTARY_CHARGE));
        CHECK(isExactDoubleEqual(atm.getMass(), 4*Core::AMU_TO_KG));
        CHECK(isExactDoubleEqual(atm.getPartCharge(), 2.0*Core::ELEMENTARY_CHARGE));
        CHECK(atm.getType() == CollisionModel::Atom::AtomType::H);
        CHECK(isExactDoubleEqual(atm.getSigma(), 0.0));
        CHECK(isExactDoubleEqual(atm.getEpsilon(), 0.0));
    }

    SECTION("Constructor with position, mass, (partial) charge and LJ parameter"){
        CollisionModel::Atom atm = CollisionModel::Atom(Core::Vector(0.5, 0.1, -0.2), 4, -1.0, 2.0, 
                                    CollisionModel::Atom::AtomType::Ar, 1.6, 0.02);

        CHECK(atm.getRelativePosition() == Core::Vector(0.5, 0.1, -0.2));
        CHECK(isExactDoubleEqual(atm.getCharge(), -1.0*Core::ELEMENTARY_CHARGE));
        CHECK(isExactDoubleEqual(atm.getMass(), 4*Core::AMU_TO_KG));
        CHECK(isExactDoubleEqual(atm.getPartCharge(), 2.0*Core::ELEMENTARY_CHARGE));
        CHECK(atm.getType() == CollisionModel::Atom::AtomType::Ar);
        CHECK(isExactDoubleEqual(atm.getSigma(), 1.6));
        CHECK(isExactDoubleEqual(atm.getEpsilon(), 0.02));
    }
}

TEST_CASE("Basic Atom setter tests", "[CollisionModels][Atom]") {
    CollisionModel::Atom atm = CollisionModel::Atom(Core::Vector(0.5, 0.1, -0.2), 4.003, -1.0, 0.2, 
                                    CollisionModel::Atom::AtomType::Ar, 1.6, 0.02);
    
    atm.setRelativePosition(Core::Vector(0.1, 0.1, 0.1));
    CHECK(atm.getRelativePosition() == Core::Vector(0.1, 0.1, 0.1));

    atm.setCharge(10);
    CHECK(isExactDoubleEqual(atm.getCharge(), 10*Core::ELEMENTARY_CHARGE));

    atm.setPartCharge(5);
    CHECK(isExactDoubleEqual(atm.getPartCharge(), 5*Core::ELEMENTARY_CHARGE));

    atm.setMass(20);
    CHECK(isExactDoubleEqual(atm.getMass(), 20*Core::AMU_TO_KG));

    atm.setType(CollisionModel::Atom::AtomType::N);
    CHECK(atm.getType() == CollisionModel::Atom::AtomType::N);

    atm.setSigma(2);
    CHECK(isExactDoubleEqual(atm.getSigma(), 2));

    atm.setEpsilon(0.4);
    CHECK(isExactDoubleEqual(atm.getEpsilon(), 0.4));

}

TEST_CASE("Atom LJ parameter calculation", "[CollisionModels][Atom]") {
    CollisionModel::Atom atm1 = CollisionModel::Atom();
    atm1.setEpsilon(0.02);
    atm1.setSigma(2.556);
    CollisionModel::Atom atm2 = CollisionModel::Atom();
    atm2.setEpsilon(0.2339);
    atm2.setSigma(3.401);

    CHECK(CollisionModel::Atom::calcLJEps(atm1, atm2) == Approx(0.0684).margin(0.0001));
    CHECK(CollisionModel::Atom::calcLJSig(atm1, atm2) == Approx(2.9785).margin(0.0001));

}

TEST_CASE("Rotation of atoms", "[CollisionModels][Atom]") {
    CollisionModel::Atom atm1 = CollisionModel::Atom();
    atm1.setRelativePosition(Core::Vector(0, 370/2, 0));
    CollisionModel::Atom atm2 = CollisionModel::Atom();
    atm2.setRelativePosition(Core::Vector(0, -370/2, 0));


    SECTION("Rotation around z"){
        atm1.rotate(Core::Vector(0, 0, M_PI/2));
        atm2.rotate(Core::Vector(0, 0, M_PI/2));

        CHECK(atm1.getRelativePosition().x() == Approx(185.0).margin(1E-13));
        CHECK(atm1.getRelativePosition().y() == Approx(0.0).margin(1E-13));
        CHECK(atm1.getRelativePosition().z() == Approx(0.0).margin(1E-13));
        CHECK(atm2.getRelativePosition().x() == Approx(-185.0).margin(1E-13));
        CHECK(atm2.getRelativePosition().y() == Approx(0.0).margin(1E-13));
        CHECK(atm2.getRelativePosition().z() == Approx(0.0).margin(1E-13));
    }

    SECTION("Rotation about y"){
        atm1.rotate(Core::Vector(0, M_PI/2, 0));
        atm2.rotate(Core::Vector(0, M_PI/2, 0));

        CHECK(atm1.getRelativePosition().x() == Approx(0.0).margin(1E-13));
        CHECK(atm1.getRelativePosition().y() == Approx(185.0).margin(1E-13));
        CHECK(atm1.getRelativePosition().z() == Approx(0.0).margin(1E-13));
        CHECK(atm2.getRelativePosition().x() == Approx(0.0).margin(1E-13));
        CHECK(atm2.getRelativePosition().y() == Approx(-185.0).margin(1E-13));
        CHECK(atm2.getRelativePosition().z() == Approx(0.0).margin(1E-13));
    }

    SECTION("Rotation about x"){
        atm1.rotate(Core::Vector(M_PI/2, 0, 0));
        atm2.rotate(Core::Vector(M_PI/2, 0, 0));

        CHECK(atm1.getRelativePosition().x() == Approx(0.0).margin(1E-13));
        CHECK(atm1.getRelativePosition().y() == Approx(0.0).margin(1E-13));
        CHECK(atm1.getRelativePosition().z() == Approx(-185.0).margin(1E-13));
        CHECK(atm2.getRelativePosition().x() == Approx(0.0).margin(1E-13));
        CHECK(atm2.getRelativePosition().y() == Approx(0.0).margin(1E-13));
        CHECK(atm2.getRelativePosition().z() == Approx(185.0).margin(1E-13));
    }

    SECTION("All three rotations"){
        atm1.rotate(Core::Vector(M_PI/2, M_PI/2, M_PI/2));
        atm2.rotate(Core::Vector(M_PI/2, M_PI/2, M_PI/2));

        CHECK(atm1.getRelativePosition().x() == Approx(0.0).margin(1E-13));
        CHECK(atm1.getRelativePosition().y() == Approx(-185.0).margin(1E-13));
        CHECK(atm1.getRelativePosition().z() == Approx(0.0).margin(1E-13));
        CHECK(atm2.getRelativePosition().x() == Approx(0.0).margin(1E-13));
        CHECK(atm2.getRelativePosition().y() == Approx(185.0).margin(1E-13));
        CHECK(atm2.getRelativePosition().z() == Approx(0.0).margin(1E-13));
    }

}
