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
 CollisionModel_Molecule.hpp

 Class to create molecules for MD Interactions Collision Model.


 ****************************/

#ifndef IDSIMF_COLLISIONMODEL_MOLECULE_H
#define IDSIMF_COLLISIONMODEL_MOLECULE_H

#include "CollisionModel_Atom.hpp"
#include "Core_constants.hpp"
#include "Core_vector.hpp"
#include "CollisionModel_MolecularStructure.hpp"
#include <vector>

namespace CollisionModel{

    class Molecule {

    public:

        // Constructors
        Molecule() = default;
        ~Molecule() = default;
        Molecule(const Core::Vector &comPos, const Core::Vector &comVel);
        Molecule(const Core::Vector &comPos, const Core::Vector &comVel, 
                 const Core::Vector &agls, std::vector<std::shared_ptr<CollisionModel::Atom>> &atms,
                 double diam);
        Molecule(const Core::Vector &comPos, const Core::Vector &comVel, std::shared_ptr<CollisionModel::MolecularStructure> structure);

        // Setter
        void setComPos(Core::Vector comPos);
        void setComVel(Core::Vector comVel);
        void setAngles(Core::Vector agls);
        void setDiameter(double diam);
        void setMolecularStructureName(std::string name);

        // Getter 
        Core::Vector& getComPos();
        Core::Vector& getComVel();
        Core::Vector& getAngles();
        bool getIsDipole() const;
        bool getIsIon() const;
        double getMass() const;
        Core::Vector getDipole() const;
        double getDipoleMag() const;
        std::size_t getAtomCount() const;
        std::vector<std::shared_ptr<CollisionModel::Atom>>& getAtoms();
        double getDiameter() const;
        std::string getMolecularStructureName() const; 

        // Member functions
        void addAtom(std::shared_ptr<CollisionModel::Atom> atm);
        void removeAtom(std::shared_ptr<CollisionModel::Atom> atm);
        void rotateMolecule();

    private:

        // Helper functions
        void calcMass();
        void calcDipole();
        void setIsDipole();
        void setIsIon();
        

        // Attributes
        Core::Vector centerOfMassPos = {0.0, 0.0, 0.0}; // Center-of-mass position of the molecule [m]
        Core::Vector centerOfMassVel = {0.0, 0.0, 0.0}; // Center-of-mass velocity of the molecule [m/s]
        Core::Vector angles = {0.0, 0.0, 0.0}; // Angles describing the xyz-rotation of the molecule w.r.t its starting configuration [rad]
        bool isDipole = false; // Flag to denote the molecule is a neutral dipole 
        bool isIon = false; // Flag to denote the molecule is an ion 
        double mass = 0.0; // Mass of the molecule [kg]
        Core::Vector dipole = {0.0, 0.0, 0.0}; // Dipole vector of the molecule [C*m]
        double dipoleMag = 0.0; // Magnitude of dipole [C]
        std::size_t atomCount = 0; // Number of atoms belonging to the molecule
        std::vector<std::shared_ptr<CollisionModel::Atom>> atoms; // Vector of all atoms belonging to this molecule 
        double diameter = 0.0; // Diameter of the molecule for collision probability [m]
        std::string molecularStructureName = "";

        
    };

}

#endif //IDSIMF_COLLISIONMODEL_MOLECULE_H
