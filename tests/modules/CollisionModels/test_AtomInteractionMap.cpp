/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2025 - Physical and Theoretical Chemistry /
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
 test_AtomInteractionMap.cpp

 Tests for the interaction map, a mapping of pairs of two atoms to other values

 ****************************/



#include "CollisionModel_AtomInteractionMap.hpp"
#include "catch.hpp"

TEST_CASE( "Test atom interaction map", "[CollisionModels][AtomInteractionMap]") {

    SECTION("AtomInteractionMap should be insertable and retrievable") {
        CollisionModel::AtomInteractionMap<double> interactionMap;

        CollisionModel::Atom atomH({1.0, 0.0, 0.0}, 1.0, 1.0);
        CollisionModel::Atom atomHe({0.0, 1.0, 0.0}, 4.0, 1.0);
        CollisionModel::Atom atomNe({0.0, 0.0, 1.0}, 20.0, 1.0);

        interactionMap.insert(atomH, atomHe, 2.0);
        CHECK(interactionMap.get(atomH, atomHe) == Approx(2.0));
        CHECK(interactionMap.get(atomHe, atomH) == Approx(2.0));

        interactionMap.insert(atomH, atomHe, 3.0);
        // currently: Re-Assignments are silently ignored
        CHECK(interactionMap.get(atomH, atomHe) == Approx(2.0));
        CHECK(interactionMap.get(atomHe, atomH) == Approx(2.0));

        interactionMap.insert(atomNe, atomH, 4.0);
        CHECK(interactionMap.get(atomH, atomNe) == Approx(4.0));
        CHECK(interactionMap.get(atomNe, atomH) == Approx(4.0));

        interactionMap.insert(atomH, atomNe, 5.0);
        CHECK(interactionMap.get(atomH, atomNe) == Approx(4.0));
        CHECK(interactionMap.get(atomNe, atomH) == Approx(4.0));

        interactionMap.insert(atomH, atomH, 10.0);
        CHECK(interactionMap.get(atomH, atomH) == Approx(10.0));
    }
}