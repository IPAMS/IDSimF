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
 test_MDInteractions.cpp

 Testing of molecular collision model with LJ-12-6 and dipole forces

 ****************************/

#include "CollisionModel_MDInteractions.hpp"
#include "CollisionModel_Molecule.hpp"
#include "CollisionModel_Atom.hpp"
#include "Core_constants.hpp"
#include "catch.hpp"
#include <iostream>

TEST_CASE("Basic test MD Interactions model", "[CollisionModels][MDInteractionsModel]") {

    double diameterHe = CollisionModel::MDInteractionsModel::DIAMETER_HE;
    CollisionModel::Molecule ion;
    CollisionModel::Molecule gas;
    double eps_Ar = 0.2339 * 4.1868 * 1E3 / Core::N_AVOGADRO; 
    double sig_Ar = 3.401 * 1E-10;
    double eps_He = 0.02 * 4.1868 * 1E3 / Core::N_AVOGADRO; 
    double sig_He = 2.556 * 1E-10;
    CollisionModel::Atom argon1 = CollisionModel::Atom(Core::Vector(0, 2.90E-10, 0), 39.948, 0, -1, 
                                    CollisionModel::Atom::AtomType::Ar, sig_Ar, eps_Ar);
    CollisionModel::Atom argon2 = CollisionModel::Atom(Core::Vector(0, -2.90E-10, 0), 39.948, 0, +1, 
                                    CollisionModel::Atom::AtomType::Ar, sig_Ar, eps_Ar);
    CollisionModel::Atom helium = CollisionModel::Atom(Core::Vector(0, 0, 0), 4.003, -1, 0, 
                                    CollisionModel::Atom::AtomType::He, sig_He, eps_He);
    ion.addAtom(&argon1);
    ion.addAtom(&argon2);
    ion.setComVel(Core::Vector(10,0,0));
    ion.setDiameter(100*4E-10);
    gas.addAtom(&helium);
    std::cout << "Vel: " << ion.getComVel().magnitude() << std::endl;

    CollisionModel::MDInteractionsModel mdSim = CollisionModel::MDInteractionsModel(1.0, 298, 4.003, 
                                                                                    diameterHe, 0.205E-30);
    //for(int i = 0; i < 200; i++)
    mdSim.modifyVelocity(ion, gas, 2e-7);

    std::cout << "Vel: " << ion.getComVel().magnitude() << std::endl;
    std::cout << "Pos: " << ion.getComPos() << std::endl;


    
}
