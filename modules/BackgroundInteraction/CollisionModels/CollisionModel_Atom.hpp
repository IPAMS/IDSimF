/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2022 - Physical and Theoretical Chemistry /
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
 CollisionModel_Atom.hpp

 Class to create individual atoms for MD Interactions Collision Model.
 Atoms are defined as a part of a parent molecule. 


 ****************************/

#ifndef IDSIMF_COLLISIONMODEL_ATOM_H
#define IDSIMF_COLLISIONMODEL_ATOM_H

#include "Core_constants.hpp"
#include "Core_vector.hpp"
#include <exception>

namespace CollisionModel{

    class Atom {

    public:

        enum class AtomType : int {C, O, N, H, He, Ar, Cl, Li, Sn, COM}; // NOTE: This might not be necessary when directly saving the LJ params 
        static AtomType from_string(std::string str); // Returns enum based on string 

        // Constructors
        Atom() = default;
        ~Atom() = default;
        Atom(const Core::Vector &relPos, double massAMU, double chargeElemCharges);
        Atom(const Core::Vector &relPos, double massAMU, double chargeElemCharges, double partChargeElemCharges);
        Atom(const Core::Vector &relPos, double massAMU, double chargeElemCharges, double partChargeElemCharges, AtomType element, double sig, double eps);

        //Setter 
        void setRelativePosition(Core::Vector relPos);
        void setMass(double massAMU);
        void setType(AtomType element);
        void setSigma(double sig);
        void setEpsilon(double eps);
        void setCharge(double chargeElemCharges);
        void setPartCharge(double partChargeElemCharges);

        //Getter
        Core::Vector& getRelativePosition();
        double getMass() const;
        AtomType getType() const;
        double getSigma() const;
        double getEpsilon() const;
        double getCharge() const;
        double getPartCharge() const;

        // Member function

        void rotate(const Core::Vector &angles); 

        //Helper functions

        static double calcLJEps(const Atom &atm1, const Atom &atm2);
        static double calcLJSig(const Atom &atm1, const Atom &atm2);
        static void rotate2D(double angle, Core::Vector& relPos); 

    private:

        // Attributes 
        Core::Vector relativePosition = Core::Vector(0.0, 0.0, 0.0); // Relative position to the center-of-mass of the parent molecule [m]
        double mass = 0.0; // Mass of the atom [kg]
        double charge = 0.0; // Charge [C]
        double partialCharge = 0.0; // Partial charge for dipole moments [C]
        AtomType type = AtomType::H; // Element of the atom 
        double sigma = 0.0; // Lennard Jones parameter sigma w.r.t its own type [m]
        double epsilon = 0.0; // Lennard Jones parameter epsilon w.r.t its own type [J]


    };

}

#endif //IDSIMF_COLLISIONMODEL_ATOM_H
