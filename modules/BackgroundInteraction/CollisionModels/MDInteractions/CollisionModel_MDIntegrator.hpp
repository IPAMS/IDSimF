/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2026 - Physical and Theoretical Chemistry /
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
 MDIntegrator.hpp

 Description

 ****************************/
#ifndef IDSIMF_MDINTEGRATOR_HPP
#define IDSIMF_MDINTEGRATOR_HPP

#include "CollisionModel_Molecule.hpp"
#include "CollisionModel_AbstractMDForceField.hpp"
#include "CollisionModel_MDTrajectoryWriter.hpp"
#include "CollisionModel_HDF5MDTrajectoryWriter.hpp"
#include "Core_randomGenerators.hpp"
#include "AppUtils_logging.hpp"

namespace CollisionModel{
    class MDIntegrator {

    public:
        enum MDIntegratorType{RK4_ADAPTIVE, RK4, LEAPFROG};

        MDIntegrator(
            std::string collisionMolecule,
            bool rotationActive,
            std::unique_ptr<AbstractMDForceField> forceField);


        void setLegacyTrajectoryWriter(const std::string& trajectoryFileName,
                         double trajectoryDistance,
                         double minimalSampleInterval,
                         unsigned int startTimeStep=0);

        void setHDF5TrajectoryWriter(const std::string& trajectoryFileName,
                                double trajectoryDistance,
                                unsigned int startTimeStep=0);

        bool leapfrogIntern(std::vector<Molecule*> moleculesPtr, double dt, double finalTime, double requiredRad);
        bool rk4Intern(std::vector<Molecule*> moleculesPtr, double dt, double finalTime, double requiredRad);
        bool rk4InternAdaptiveStep(std::vector<Molecule*> moleculesPtr, double dt, double finalTime, int maximumTimeSteps, double requiredRad, double tolerance, bool ionIsFrozen=false);

    protected:

        struct {
            bool modelRecordsTrajectory = false;
            bool recordingActive = false;
            double trajectoryDistance = 0.0;
            unsigned int recordTrajectoryStartTimeStep = 0;
        } legacyTWriterConf_, hdf5TWriterConf_;

        double initRotation(Molecule& mole, double temperature_K);
        double calcRotEnergy(Core::Vector omega, Core::Matrix3 I);

        std::string collisionMolecule_ = "";
        bool rotationActive_ = false;

        std::unique_ptr<MDTrajectoryWriter> legacyTrajectoryWriter_ = nullptr;
        std::unique_ptr<HDF5MDTrajectoryWriter> hdf5TrajectoryWriter_ = nullptr;

        std::unique_ptr<AbstractMDForceField> forceField_; ///< The molecular force field to use
    };
}

#endif //IDSIMF_MDINTEGRATOR_HPP