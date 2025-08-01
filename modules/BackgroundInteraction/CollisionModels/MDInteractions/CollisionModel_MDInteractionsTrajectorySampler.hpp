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
 CollisionModel_MDInteractionsTrajectorySampler.hpp

 Explicit Molecular Dynamics Collision Trajectory Sampler / Collision Trajectory Calculation

 Allows to calculate explicit MD collision trajectories with defined initial and boundary conditions primarily for
 detailed analysis of molecular collision dynamics
 ****************************/

#ifndef IDSIMF_COLLISIONMODEL_MDINTERACTIONS_TRAJECTORY_SAMPLER_H
#define IDSIMF_COLLISIONMODEL_MDINTERACTIONS_TRAJECTORY_SAMPLER_H

#include "Core_constants.hpp"
#include "CollisionModel_MDTrajectoryWriter.hpp"
#include "CollisionModel_AbstractMDForceField.hpp"
#include "CollisionModel_Molecule.hpp"
#include "RS_AbstractReaction.hpp"
#include "AppUtils_logging.hpp"
#include <string>

namespace CollisionModel{

    enum MDIntegratorType{RK4_ADAPTIVE, RK4, LEAPFROG};

    class MDInteractionsTrajectorySampler {

    public:
        constexpr static double DIAMETER_N2 = 3.64e-10;
        constexpr static double DIAMETER_HE = 2.80e-10;

        MDInteractionsTrajectorySampler() = default;

        MDInteractionsTrajectorySampler(
            std::unique_ptr<AbstractMDForceField> forceField,
            std::unordered_map<std::string, std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection,
            AppUtils::logger_ptr logger);

        void setTrajectoryWriter(const std::string& trajectoryFileName,
                                 double minimalSampleInterval=0);

        void calculateTrajectory(Core::Particle& particle,
                                 Core::Vector particleRotationAngles,
                                 std::string collisionMolecule,
                                 Core::Vector collisionParticleStartPosition,
                                 Core::Vector collisionParticleStartVelocity,
                                 Core::Vector collisionParticleRotationAngles,
                                 double integrationTime,
                                 double subTimeStep,
                                 int maximumSteps,
                                 bool ionIsFrozen,
                                 MDIntegratorType integratorType = RK4_ADAPTIVE);

        bool leapfrogIntern(std::vector<CollisionModel::Molecule*> moleculesPtr, double dt, double finalTime, int maximumSteps, bool ionIsFrozen);
        bool rk4Intern(std::vector<CollisionModel::Molecule*> moleculesPtr, double dt, double finalTime, int maximumSteps, bool ionIsFrozen);
        bool rk4InternAdaptiveStep(std::vector<CollisionModel::Molecule*> moleculesPtr, double dt, double finalTime, int maximumSteps, bool ionIsFrozen);

    private:

        std::string currentCollisionMolecule_ = "";

        std::unique_ptr<MDTrajectoryWriter> trajectoryWriter_ = nullptr;
        std::unique_ptr<AbstractMDForceField> forceField_; ///< The molecular force field to use
        std::unordered_map<std::string,  std::shared_ptr<MolecularStructure>> molecularStructureCollection_; ///< collection of all available molecular structures
        AppUtils::logger_ptr logger_ = nullptr;
    };
}

#endif //IDSIMF_COLLISIONMODEL_MDINTERACTIONS_TRAJECTORY_SAMPLER_H
