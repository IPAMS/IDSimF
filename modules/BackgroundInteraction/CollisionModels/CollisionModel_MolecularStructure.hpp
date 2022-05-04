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
 CollisionModel_MolecularStructure.hpp

 Class to create molecules for MD Interactions Collision Model.


 ****************************/

#ifndef IDSIMF_COLLISIONMODEL_MOLECULARSTRUCTURE_H
#define IDSIMF_COLLISIONMODEL_MOLECULARSTRUCTURE_H

#include "CollisionModel_Atom.hpp"
#include "Core_constants.hpp"
#include "Core_vector.hpp"
#include <vector>
#include <unordered_map>
#include <string>
#include <memory>

namespace CollisionModel{

    class MolecularStructure {

    public:
        // Hash map to collect molecular structures so they can be constructed in the collision model later 
        static std::unordered_map<std::string,  std::unique_ptr<MolecularStructure>> molecularStructureCollection;

        // Constructors
        MolecularStructure() = default;
        ~MolecularStructure() = default;
        MolecularStructure(std::vector<CollisionModel::Atom*> atms, double diam);

        // Setter
        void setDiameter(double diam);

        // Getter
        bool getIsDipole() const;
        bool getIsIon() const;
        double getMass() const;
        Core::Vector getDipole() const;
        double getDipoleMag() const;
        std::size_t getAtomCount() const;
        std::vector<CollisionModel::Atom*> getAtoms() const;
        double getDiameter() const;

        // Member functions
        void addAtom(CollisionModel::Atom* atm);
        void removeAtom(CollisionModel::Atom* atm);

        static std::unordered_map<std::string, std::unique_ptr<MolecularStructure>> createCollection();
        

    private:

        // Helper functions
        void calcMass();
        void calcDipole();
        void setIsDipole();
        void setIsIon();
        

        // Attributes
        bool isDipole = false; // Flag to denote the molecule is a neutral dipole 
        bool isIon = false; // Flag to denote the molecule is an ion 
        double mass = 0.0; // Mass of the molecule [kg]
        Core::Vector dipole = {0.0, 0.0, 0.0}; // Dipole vector of the molecule [C*m]
        double dipoleMag = 0.0; // Magnitude of dipole [C]
        std::size_t atomCount = 0; // Number of atoms belonging to the molecule
        std::vector<CollisionModel::Atom*> atoms = {}; // Vector of all atoms belonging to this molecule 
        double diameter = 0.0; // Diameter of the molecule for collision probability [m]

        
    };

}

#endif //IDSIMF_COLLISIONMODEL_MOLECULARSTRUCTURE_H
