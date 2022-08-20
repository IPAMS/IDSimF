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
 test_Molecule.cpp

 Testing of Molecule class for use in MD Interactions Collision Model.

 ****************************/

#include "CollisionModel_Molecule.hpp"
#include "Core_constants.hpp"
#include "catch.hpp"
#include "test_util.hpp"
#include "CollisionModel_MolecularStructure.hpp"
#include <iostream>
#include <memory>

TEST_CASE("Basic test Molecule creation", "[CollisionModels][Molecule]") {

    SECTION("Default constructor"){
        CollisionModel::Molecule mole = CollisionModel::Molecule();

        CHECK(mole.getComPos() == Core::Vector(0.0, 0.0, 0.0));
        CHECK(mole.getComVel() == Core::Vector(0.0, 0.0, 0.0));
        CHECK(mole.getAngles() == Core::Vector(0.0, 0.0, 0.0));
        CHECK(mole.getDipole() == Core::Vector(0.0, 0.0, 0.0));
        CHECK(mole.getIsDipole() == false);
        CHECK(mole.getIsIon() == false);
        CHECK(isExactDoubleEqual(mole.getMass(), 0.0*Core::AMU_TO_KG));
        CHECK(isExactDoubleEqual(mole.getDipoleMag(), 0.0));
        CHECK(mole.getAtomCount() == 0);
        CHECK(mole.getAtoms().empty() == true);
    }

    SECTION("Constructor with position and velocity"){
        CollisionModel::Molecule mole = CollisionModel::Molecule(Core::Vector(0.5, 1.0, -0.3), 
                                                                    Core::Vector(-1.5, 0.3, 0.33));

        CHECK(mole.getComPos() == Core::Vector(0.5, 1.0, -0.3));
        CHECK(mole.getComVel() == Core::Vector(-1.5, 0.3, 0.33));
        CHECK(mole.getAngles() == Core::Vector(0.0, 0.0, 0.0));
        CHECK(mole.getDipole() == Core::Vector(0.0, 0.0, 0.0));
        CHECK(mole.getIsDipole() == false);
        CHECK(mole.getIsIon() == false);
        CHECK(isExactDoubleEqual(mole.getMass(), 0.0*Core::AMU_TO_KG));
        CHECK(isExactDoubleEqual(mole.getDipoleMag(), 0.0));
        CHECK(mole.getAtomCount() == 0);
        CHECK(mole.getAtoms().empty() == true);
    }

    SECTION("Constructor with position, velocity and atoms"){
        CollisionModel::Atom atm1 = CollisionModel::Atom();
        atm1.setRelativePosition(Core::Vector(0, 370/2, 0));
        atm1.setMass(1);
        CollisionModel::Atom atm2 = CollisionModel::Atom();
        atm2.setRelativePosition(Core::Vector(0, -370/2, 0));
        atm2.setMass(3);

        std::vector<std::shared_ptr<CollisionModel::Atom>> atoms = {
            std::make_shared<CollisionModel::Atom>(std::move(atm1)), std::make_shared<CollisionModel::Atom>(std::move(atm2))};
        CollisionModel::Molecule mole = CollisionModel::Molecule(Core::Vector(0.5, 1.0, -0.3), 
                                                                    Core::Vector(-1.5, 0.3, 0.33),
                                                                    Core::Vector(0.0, 0.0, 0.0),
                                                                    atoms, 0.5);

        CHECK(mole.getComPos() == Core::Vector(0.5, 1.0, -0.3));
        CHECK(mole.getComVel() == Core::Vector(-1.5, 0.3, 0.33));
        CHECK(mole.getAngles() == Core::Vector(0.0, 0.0, 0.0));
        CHECK(mole.getDipole() == Core::Vector(0.0, 0.0, 0.0));
        CHECK(mole.getIsDipole() == false);
        CHECK(mole.getIsIon() == false);
        CHECK(isExactDoubleEqual(mole.getMass(), 4*Core::AMU_TO_KG));
        CHECK(isExactDoubleEqual(mole.getDipoleMag(), 0.0));
        CHECK(isExactDoubleEqual(mole.getAtomCount(), 2));
        CHECK(mole.getAtoms().empty() == false);
        CHECK(isExactDoubleEqual(mole.getDiameter(), 0.5));
    }

    SECTION("Constructor with position, velocity and atoms (dipole and ion)"){
        CollisionModel::Atom atm1 = CollisionModel::Atom();
        atm1.setRelativePosition(Core::Vector(0, 370/2, 0));
        atm1.setMass(1);
        atm1.setPartCharge(-0.5);
        CollisionModel::Atom atm2 = CollisionModel::Atom();
        atm2.setRelativePosition(Core::Vector(0, -370/2, 0));
        atm2.setMass(3);
        atm2.setCharge(1);
        atm2.setPartCharge(0.5);

        
        std::vector<std::shared_ptr<CollisionModel::Atom>> atoms = {
            std::make_shared<CollisionModel::Atom>(std::move(atm1)), std::make_shared<CollisionModel::Atom>(std::move(atm2))};
        CollisionModel::Molecule mole = CollisionModel::Molecule(Core::Vector(0.5, 1.0, -0.3), 
                                                                    Core::Vector(-1.5, 0.3, 0.33),
                                                                    Core::Vector(0.0, 0.0, 0.0),
                                                                    atoms, 0.5);

        CHECK(mole.getComPos() == Core::Vector(0.5, 1.0, -0.3));
        CHECK(mole.getComVel() == Core::Vector(-1.5, 0.3, 0.33));
        CHECK(mole.getAngles() == Core::Vector(0.0, 0.0, 0.0));
        CHECK(mole.getDipole().y() == Approx(-2.96E-17).margin(0.001E-15));
        CHECK(mole.getIsDipole() == true);
        CHECK(mole.getIsIon() == true);
        CHECK(isExactDoubleEqual(mole.getMass(), 4*Core::AMU_TO_KG));
        CHECK(mole.getDipoleMag() == Approx(2.96E-17).margin(0.001E-15));
        CHECK(isExactDoubleEqual(mole.getAtomCount(), 2));
        CHECK(mole.getAtoms().empty() == false);
        CHECK(isExactDoubleEqual(mole.getDiameter(), 0.5));
    }

    SECTION("Constructor with position, velocity and atoms (non-dipole and ion)"){
        CollisionModel::Atom atm1 = CollisionModel::Atom();
        atm1.setRelativePosition(Core::Vector(0, 370/2, 0));
        atm1.setMass(1);
        CollisionModel::Atom atm2 = CollisionModel::Atom();
        atm2.setRelativePosition(Core::Vector(0, -370/2, 0));
        atm2.setMass(3);
        atm2.setCharge(1);

        std::vector<std::shared_ptr<CollisionModel::Atom>> atoms = {
            std::make_shared<CollisionModel::Atom>(std::move(atm1)), std::make_shared<CollisionModel::Atom>(std::move(atm2))};
        CollisionModel::Molecule mole = CollisionModel::Molecule(Core::Vector(0.5, 1.0, -0.3), 
                                                                    Core::Vector(-1.5, 0.3, 0.33),
                                                                    Core::Vector(0.0, 0.0, 0.0),
                                                                    atoms, 0.5);

        CHECK(mole.getComPos() == Core::Vector(0.5, 1.0, -0.3));
        CHECK(mole.getComVel() == Core::Vector(-1.5, 0.3, 0.33));
        CHECK(mole.getAngles() == Core::Vector(0.0, 0.0, 0.0));
        CHECK(mole.getDipole().y() == Approx(0.0).margin(1E-13));
        CHECK(mole.getIsDipole() == false);
        CHECK(mole.getIsIon() == true);
        CHECK(isExactDoubleEqual(mole.getMass(), 4*Core::AMU_TO_KG));
        CHECK(mole.getDipoleMag() == Approx(0.0).margin(1E-13));
        CHECK(isExactDoubleEqual(mole.getAtomCount(), 2));
        CHECK(mole.getAtoms().empty() == false);
        CHECK(isExactDoubleEqual(mole.getDiameter(), 0.5));
    }

    SECTION("Constructor with position, velocity and atoms (non-dipole and non-ion)"){
        CollisionModel::Atom atm1 = CollisionModel::Atom();
        atm1.setRelativePosition(Core::Vector(0, 370/2, 0));
        atm1.setMass(1);
        CollisionModel::Atom atm2 = CollisionModel::Atom();
        atm2.setRelativePosition(Core::Vector(0, -370/2, 0));
        atm2.setMass(3);

        std::vector<std::shared_ptr<CollisionModel::Atom>> atoms = {
            std::make_shared<CollisionModel::Atom>(std::move(atm1)), std::make_shared<CollisionModel::Atom>(std::move(atm2))};
        CollisionModel::Molecule mole = CollisionModel::Molecule(Core::Vector(0.5, 1.0, -0.3), 
                                                                    Core::Vector(-1.5, 0.3, 0.33),
                                                                    Core::Vector(0.0, 0.0, 0.0),
                                                                    atoms, 0.5);

        CHECK(mole.getComPos() == Core::Vector(0.5, 1.0, -0.3));
        CHECK(mole.getComVel() == Core::Vector(-1.5, 0.3, 0.33));
        CHECK(mole.getAngles() == Core::Vector(0.0, 0.0, 0.0));
        CHECK(mole.getDipole().y() == Approx(0.0).margin(1E-13));
        CHECK(mole.getIsDipole() == false);
        CHECK(mole.getIsIon() == false);
        CHECK(isExactDoubleEqual(mole.getMass(), 4*Core::AMU_TO_KG));
        CHECK(mole.getDipoleMag() == Approx(0.0).margin(1E-13));
        CHECK(isExactDoubleEqual(mole.getAtomCount(), 2));
        CHECK(mole.getAtoms().empty() == false);
        CHECK(isExactDoubleEqual(mole.getDiameter(), 0.5));
    }

    SECTION("Constructor through molecular structure object"){
        CollisionModel::Atom atm1 = CollisionModel::Atom();
        atm1.setRelativePosition(Core::Vector(0, 370/2, 0));
        atm1.setMass(1);
        CollisionModel::Atom atm2 = CollisionModel::Atom();
        atm2.setRelativePosition(Core::Vector(0, -370/2, 0));
        atm2.setMass(3);
        std::vector<std::shared_ptr<CollisionModel::Atom>> atoms = {
            std::make_shared<CollisionModel::Atom>(std::move(atm1)), std::make_shared<CollisionModel::Atom>(std::move(atm2))};

        std::shared_ptr<CollisionModel::MolecularStructure> molstr = std::make_shared<CollisionModel::MolecularStructure>(atoms, 0.5, "Test");

        CollisionModel::Molecule mole = CollisionModel::Molecule(Core::Vector(0.5, 1.0, -0.3), 
                                                                    Core::Vector(-1.5, 0.3, 0.33),
                                                                    molstr);

        CHECK(mole.getComPos() == Core::Vector(0.5, 1.0, -0.3));
        CHECK(mole.getComVel() == Core::Vector(-1.5, 0.3, 0.33));
        CHECK(mole.getAngles() == Core::Vector(0.0, 0.0, 0.0));
        CHECK(mole.getDipole().y() == Approx(0.0).margin(1E-13));
        CHECK(mole.getIsDipole() == false);
        CHECK(mole.getIsIon() == false);
        CHECK(isExactDoubleEqual(mole.getMass(), 4*Core::AMU_TO_KG));
        CHECK(mole.getDipoleMag() == Approx(0.0).margin(1E-13));
        CHECK(isExactDoubleEqual(mole.getAtomCount(), 2));
        CHECK(mole.getAtoms().empty() == false);
        CHECK(isExactDoubleEqual(mole.getDiameter(), 0.5));
    }

    SECTION("Constructor through molecular structure object, leaving structure untouched"){
        CollisionModel::Atom atm1 = CollisionModel::Atom();
        atm1.setRelativePosition(Core::Vector(0, 370/2, 0));
        atm1.setMass(1);
        CollisionModel::Atom atm2 = CollisionModel::Atom();
        atm2.setRelativePosition(Core::Vector(0, -370/2, 0));
        atm2.setMass(3);
        std::vector<std::shared_ptr<CollisionModel::Atom>> atoms = {
            std::make_shared<CollisionModel::Atom>(std::move(atm1)), std::make_shared<CollisionModel::Atom>(std::move(atm2))};

        std::shared_ptr<CollisionModel::MolecularStructure> molstr = std::make_shared<CollisionModel::MolecularStructure>(atoms, 0.5, "Test");

        CollisionModel::Molecule mole = CollisionModel::Molecule(Core::Vector(0.5, 1.0, -0.3), 
                                                                    Core::Vector(-1.5, 0.3, 0.33),
                                                                    molstr);

        CHECK(mole.getComPos() == Core::Vector(0.5, 1.0, -0.3));
        CHECK(mole.getComVel() == Core::Vector(-1.5, 0.3, 0.33));
        CHECK(mole.getAngles() == Core::Vector(0.0, 0.0, 0.0));
        CHECK(mole.getDipole().y() == Approx(0.0).margin(1E-13));
        CHECK(mole.getIsDipole() == false);
        CHECK(mole.getIsIon() == false);
        CHECK(isExactDoubleEqual(mole.getMass(), 4*Core::AMU_TO_KG));
        CHECK(mole.getDipoleMag() == Approx(0.0).margin(1E-13));
        CHECK(isExactDoubleEqual(mole.getAtomCount(), 2));
        CHECK(mole.getAtoms().empty() == false);
        CHECK(isExactDoubleEqual(mole.getDiameter(), 0.5));

        mole.getAtoms().at(0)->getRelativePosition().y(2);
        CHECK(molstr->getAtoms().at(0)->getRelativePosition().y() == Approx(370/2).margin(1E-13));
        CHECK(mole.getAtoms().at(0)->getRelativePosition().y() == Approx(2).margin(1E-13));
    }
}


TEST_CASE("Basic Molecule setter tests", "[CollisionModels][Molecule]") {
    CollisionModel::Atom atm1 = CollisionModel::Atom();
    CollisionModel::Atom atm2 = CollisionModel::Atom();

        std::vector<std::shared_ptr<CollisionModel::Atom>> atoms = {
            std::make_shared<CollisionModel::Atom>(std::move(atm1)), std::make_shared<CollisionModel::Atom>(std::move(atm2))};
    CollisionModel::Molecule mole = CollisionModel::Molecule(Core::Vector(0.5, 1.0, -0.3), 
                                                                Core::Vector(-1.5, 0.3, 0.33),
                                                                Core::Vector(0.0, 0.0, 0.0),
                                                                atoms, 0.5);
    
    mole.setComPos(Core::Vector(0.1, 0.1, 0.1));
    CHECK(mole.getComPos() == Core::Vector(0.1, 0.1, 0.1));

    mole.setComVel(Core::Vector(0.1, 0.1, 0.1));
    CHECK(mole.getComVel() == Core::Vector(0.1, 0.1, 0.1));

    mole.setAngles(Core::Vector(0.1, 0.1, 0.1));
    CHECK(mole.getAngles() == Core::Vector(0.1, 0.1, 0.1));

    mole.setDiameter(0.2);
    CHECK(isExactDoubleEqual(mole.getDiameter(), 0.2));
    

}

TEST_CASE("Molecule add/remove atoms tests", "[CollisionModels][Molecule]") {
    CollisionModel::Atom atm1 = CollisionModel::Atom();
    CollisionModel::Atom atm2 = CollisionModel::Atom();

    std::vector<std::shared_ptr<CollisionModel::Atom>> atoms = {};
    CollisionModel::Molecule mole = CollisionModel::Molecule(Core::Vector(0.5, 1.0, -0.3), 
                                                                Core::Vector(-1.5, 0.3, 0.33),
                                                                Core::Vector(0.0, 0.0, 0.0),
                                                                atoms, 0.5);
                                                   
    mole.addAtom(std::make_shared<CollisionModel::Atom>(std::move(atm1)));
    CHECK(isExactDoubleEqual(mole.getAtomCount(), 1));
    mole.addAtom(std::make_shared<CollisionModel::Atom>(std::move(atm2)));
    CHECK(isExactDoubleEqual(mole.getAtomCount(), 2));
    mole.removeAtom(std::make_shared<CollisionModel::Atom>(std::move(atm2)));
    CHECK(isExactDoubleEqual(mole.getAtomCount(), 1));
    
}

TEST_CASE("Molecule ability to change properties of atoms", "[CollisionModels][Molecule]") {
    CollisionModel::Atom atm1 = CollisionModel::Atom();

    std::vector<std::shared_ptr<CollisionModel::Atom>> atoms = {};
    CollisionModel::Molecule mole = CollisionModel::Molecule(Core::Vector(0.5, 1.0, -0.3), 
                                                                Core::Vector(-1.5, 0.3, 0.33),
                                                                Core::Vector(0.0, 0.0, 0.0),
                                                                atoms, 0.5);
                                                    
    mole.addAtom(std::make_shared<CollisionModel::Atom>(std::move(atm1)));
    for(auto& atom : mole.getAtoms()){
        atom->setMass(5);
    }
    CHECK(isExactDoubleEqual(mole.getAtoms().at(0)->getMass(), 5*Core::AMU_TO_KG));
    
    
}