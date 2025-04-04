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
#include "CollisionModel_AbstractCollisionModel.hpp"
#include "CollisionModel_SpatialFieldFunctions.hpp"
#include "CollisionModel_AbstractMDForceField.hpp"
#include "CollisionModel_MathFunctions.hpp"
#include "CollisionModel_Molecule.hpp"
#include "RS_AbstractReaction.hpp"
#include "appUtils_logging.hpp"
#include <cstdio>
#include <functional>
#include <string>

namespace CollisionModel{

    class MDInteractionsTrajectorySampler {

    public:
        constexpr static double DIAMETER_N2 = 3.64e-10;
        constexpr static double DIAMETER_HE = 2.80e-10;

        MDInteractionsTrajectorySampler() = default;

        MDInteractionsTrajectorySampler(
            double collisionGasDiameterM, 
            std::string collisionMolecule,
            double integrationTime,
            double subTimeStep,
            double collisionRadiusScaling,
            //double angleThetaScaling,
            //double spawnRadius,
            std::unique_ptr<AbstractMDForceField> forceField,
            std::unordered_map<std::string, std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection, 
            Core::Vector startPosition, 
            Core::Vector startVelocity);

        void setTrajectoryWriter(const std::string& trajectoryFileName,
                                 double trajectoryDistance,
                                 unsigned int startTimeStep=0,
                                 double minimalSampleInterval=0);


        void writeTrajectory(double distance, Core::Vector positionBgMolecule, Core::Vector velocityBgMolecule, 
                        std::vector<Core::Vector> forceMolecules, bool endOfTrajectory, std::ofstream* file, double time, double dt,
                        Core::Vector positionMolecule);

        bool leapfrogIntern(std::vector<CollisionModel::Molecule*> moleculesPtr, double dt, double finalTime, double requiredRad);

        bool rk4Intern(std::vector<CollisionModel::Molecule*> moleculesPtr, double dt, double finalTime, double requiredRad);

        bool rk4InternAdaptiveStep(std::vector<CollisionModel::Molecule*> moleculesPtr, double dt, double finalTime, double requiredRad);

        void initializeModelParticleParameters(Core::Particle& ion) const;

        void updateModelParticleParameters(Core::Particle& ion) const;

        void updateModelTimestepParameters(unsigned int timestep, double time);

        void modifyAcceleration(Core::Vector& acceleration,
                                        Core::Particle& particle,
                                        double dt);
                                    
        void modifyVelocity(Core::Particle& particle);

        void modifyVelocity(Core::Particle& particle,
                                    double dt);

        void modifyPosition(Core::Vector& position,
                                    Core::Particle& particle,
                                    double dt);


    private:
        //double collisionGasMass_kg_ = 0.0;    ///< mass of the neutral colliding gas particles in kg
        double collisionGasDiameter_m_ = 0.0; ///< effective collision diameter of the neutral collision gas particles in m
        std::string collisionMolecule_ = ""; ///< particle identifier of the collision gas
        double integrationTime_ = 0.0; ///< integration time of the sub-integrator 
        double subTimeStep_ = 0.0; ///< step size of the sub-integrator 
        double collisionRadiusScaling_ = 0.0; ///< scaling for the radius of the collision sphere
        //double angleThetaScaling_ = 0.0; ///<  scaling for the angle theta
        //double spawnRadius_ = 0.0; ///< radius of the spawn sphere for the background gas particle
        double trajectoryDistance_ = 0.0; ///< distance at which the trajectory recording begins
        bool trajectoryRecordingActive_ = false; 
        bool modelRecordsTrajectories_ = false; ///< flag if trajectory will be recorded
        unsigned int recordTrajectoryStartTimeStep_ = 0; ///< time step at which the trajectory recording begins
        double recordTrajectoryMinimalSampleInterval_ = 0; ///< minimal time interval between trajectory samples
        double nextTrajectorySampleTime_ = 0; ///< next time for a written trajectory sample
        std::unique_ptr<std::ofstream> trajectoryOutputStream_;

        std::unique_ptr<AbstractMDForceField> forceField_; ///< The molecular force field to use
        //std::function<double(Core::Vector&)> pressureFunction_ = nullptr; ///< a spatial pressure function
        //std::function<Core::Vector(Core::Vector&)> velocityFunction_ = nullptr; ///< a spatial velocity function
        //std::function<double(const Core::Vector&)>temperatureFunction_ = nullptr;  ///< Spatial temperature function
        std::function<void(RS::CollisionConditions, Core::Particle&)> afterCollisionActionFunction_ = nullptr;
        ///< Function with things to do after a collision (e.g. collision based chemical reactions)
        std::unordered_map<std::string,  std::shared_ptr<MolecularStructure>> molecularStructureCollection_; ///< collection of all available molecular structures 
        Core::Vector startPosition_ = Core::Vector(0.0, 0.0, 0.0); 
        Core::Vector startVelocity_ = Core::Vector(0.0, 0.0, 0.0);

        void writeTrajectoryDelimiter_(std::ofstream* file);
    };

}

#endif //IDSIMF_COLLISIONMODEL_MDINTERACTIONS_TRAJECTORY_SAMPLER_H
