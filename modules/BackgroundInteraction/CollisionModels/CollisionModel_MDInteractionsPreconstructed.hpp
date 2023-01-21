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
 CollisionModel_MDInteractions.hpp

 Molecular collision model including LJ-12-6 potential and additionally ion-induced dipole as well as ion-permanent
 dipole interactions. Initial collisions are "constructed" through the hard sphere approach.

 This model follows the modelling of the HS1 collision model by David Manura,
 for SIMION 8.0 (Scientific Instrument Services, Inc.).
 https://simion.com/
 https://simion.com/info/collision_model_hs1.html

 Earlier hard sphere collision models:
1. Appelhans, A.D., Dahl, D.A.: Measurement of external ion injection and trapping efficiency in the ion
 trap mass spectrometer and comparison with a predictive model.
 International Journal of Mass Spectrometry. 216, 269–284 (2002). https://doi.org/10.1016/S1387-3806(02)00627-9
2. Ding, L., Sudakov, M., Kumashiro, S.: A simulation study of the digital ion trap mass spectrometer.
 International Journal of Mass Spectrometry. 221, 117–138 (2002). https://doi.org/10.1016/S1387-3806(02)00921-1

 ****************************/

#ifndef IDSIMF_COLLISIONMODEL_MDINTERACTIONS_PRECONSTRUCTED_H
#define IDSIMF_COLLISIONMODEL_MDINTERACTIONS_PRECONSTRUCTED_H

#include "Core_constants.hpp"
#include "CollisionModel_AbstractCollisionModel.hpp"
#include "CollisionModel_SpatialFieldFunctions.hpp"
#include "CollisionModel_MathFunctions.hpp"
#include "CollisionModel_Molecule.hpp"
#include "RS_AbstractReaction.hpp"
#include "appUtils_logging.hpp"
#include <cstdio>
#include <functional>
#include <string>

namespace CollisionModel{

    class MDInteractionsModelPreconstructed : public AbstractCollisionModel {

    public:
        constexpr static double DIAMETER_N2 = 3.64e-10;
        constexpr static double DIAMETER_HE = 2.80e-10;

        MDInteractionsModelPreconstructed() = default;
        
        MDInteractionsModelPreconstructed(
            double staticPressure,
            double staticTemperature,
            double collisionGasMassAmu,
            double collisionGasDiameterM, 
            double collisionGasPolarizabilityM3,
            std::string collisionMolecule,
            double integrationTime,
            double subTimeStep,
            double collisionRadiusScaling,
            double angleThetaScaling, 
            double spawnRadius,
            std::unordered_map<std::string,
            std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection);

        MDInteractionsModelPreconstructed(
            double staticPressure,
            double staticTemperature,
            double collisionGasMassAmu,
            double collisionGasDiameterM, 
            double collisionGasPolarizabilityM3,
            std::string collisionMolecule,
            double integrationTime,
            double subTimeStep,
            double collisionRadiusScaling,
            double angleThetaScaling, 
            double spawnRadius,
            std::unordered_map<std::string,
            std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection, 
            Core::Vector startPosition, 
            Core::Vector startRotation);

        MDInteractionsModelPreconstructed(
            std::function<double(Core::Vector& location)> pressureFunction,
            std::function<Core::Vector(Core::Vector& location)> velocityFunction,
            double StaticTemperature,
            double collisionGasMassAmu,
            double collisionGasDiameterM, 
            double collisionGasPolarizabilityM3,
            std::string collisionMolecule,
            double integrationTime,
            double subTimeStep,
            double collisionRadiusScaling,
            double angleThetaScaling,
            double spawnRadius,
            std::unordered_map<std::string, std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection);

        MDInteractionsModelPreconstructed(
            std::function<double(Core::Vector& location)> pressureFunction,
            std::function<Core::Vector(Core::Vector& location)> velocityFunction,
            double StaticTemperature,
            double collisionGasMassAmu,
            double collisionGasDiameterM, 
            double collisionGasPolarizabilityM3,
            std::string collisionMolecule,
            double integrationTime,
            double subTimeStep,
            double collisionRadiusScaling,
            double angleThetaScaling,
            double spawnRadius,
            std::unordered_map<std::string, std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection,
            Core::Vector startPosition, 
            Core::Vector startRotation);

        MDInteractionsModelPreconstructed(
            std::function<double(Core::Vector& location)> pressureFunction,
            std::function<Core::Vector(Core::Vector& location)> velocityFunction,
            std::function<double(const Core::Vector&)> temperatureFunction,
            double collisionGasMassAmu,
            double collisionGasDiameterM, 
            double collisionGasPolarizabilityM3,
            std::string collisionMolecule,
            double integrationTime,
            double subTimeStep,
            double collisionRadiusScaling,
            double angleThetaScaling,
            double spawnRadius,
            std::unordered_map<std::string, std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection);

        MDInteractionsModelPreconstructed(
            std::function<double(Core::Vector& location)> pressureFunction,
            std::function<Core::Vector(Core::Vector& location)> velocityFunction,
            std::function<double(const Core::Vector&)> temperatureFunction,
            double collisionGasMassAmu,
            double collisionGasDiameterM, 
            double collisionGasPolarizabilityM3,
            std::string collisionMolecule,
            double integrationTime,
            double subTimeStep,
            double collisionRadiusScaling,
            double angleThetaScaling,
            double spawnRadius,
            std::unordered_map<std::string, std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection, 
            Core::Vector startPosition, 
            Core::Vector startRotation);

        void setTrajectoryWriter(const std::string& trajectoryFileName,
                                 double trajectoryDistance,
                                 int startTimeStep=0);

        double calcSign(double value);

        Core::Vector rotate(const Core::Vector &angles, Core::Vector& position);  

        void writeTrajectory(double distance, Core::Vector positionBgMolecule, Core::Vector velocityBgMolecule, 
                        std::vector<Core::Vector> forceMolecules, bool endOfTrajectory, std::ofstream* file, double time, double dt);

        bool leapfrogIntern(std::vector<CollisionModel::Molecule*> moleculesPtr, double dt, double finalTime, double requiredRad);

        bool rk4Intern(std::vector<CollisionModel::Molecule*> moleculesPtr, double dt, double finalTime, double requiredRad);

        bool rk4InternAdaptiveStep(std::vector<CollisionModel::Molecule*> moleculesPtr, double dt, double finalTime, double requiredRad);

        void forceFieldMD(std::vector<CollisionModel::Molecule*>& moleculesPtr, std::vector<Core::Vector>& forceMolecules);

        void initializeModelParticleParameters(Core::Particle& ion) const;

        void updateModelParticleParameters(Core::Particle& ion) const;

        void updateModelTimestepParameters(int timestep, double time);

        void modifyAcceleration(Core::Vector& acceleration,
                                        Core::Particle& particle,
                                        double dt);
        void modifyVelocity(Core::Particle& particle,
                                    double dt);

        void modifyPosition(Core::Vector& position,
                                    Core::Particle& particle,
                                    double dt);


    private:
        double collisionGasMass_kg_ = 0.0;    ///< mass of the neutral colliding gas particles in kg
        double collisionGasDiameter_m_ = 0.0; ///< effective collision diameter of the neutral collision gas particles in m
        double collisionGasPolarizability_m3_ = 0.0; ///< polarizability of the collision gas in m^3
        std::string collisionMolecule_ = ""; ///< particle identifier of the collision gas
        double integrationTime_ = 0.0; ///< integration time of the sub-integrator 
        double subTimeStep_ = 0.0; ///< step size of the sub-integrator 
        double collisionRadiusScaling_ = 0.0; ///< scaling for the radius of the collision sphere 
        double angleThetaScaling_ = 0.0; ///<  scaling for the angle theta 
        double spawnRadius_ = 0.0; ///< radius of the spawn sphere for the background gas particle
        double trajectoryDistance_ = 0.0; ///< distance at which the trajectory recording begins
        bool trajectoryRecordingActive_ = false; 
        bool modelRecordsTrajectories_ = false; ///< flag if trajectory will be recorded
        int recordTrajectoryStartTimeStep_ = 0; ///< time step at which the trajectory recording begins
        std::unique_ptr<std::ofstream> trajectoryOutputStream_;

        std::function<double(Core::Vector&)> pressureFunction_ = nullptr; ///< a spatial pressure function
        std::function<Core::Vector(Core::Vector&)> velocityFunction_ = nullptr; ///< a spatial velocity function
        std::function<double(const Core::Vector&)>temperatureFunction_ = nullptr;  ///< Spatial temperature function
        std::function<void(RS::CollisionConditions, Core::Particle&)> afterCollisionActionFunction_ = nullptr;
        ///< Function with things to do after a collision (e.g. collision based chemical reactions)
        std::unordered_map<std::string,  std::shared_ptr<MolecularStructure>> molecularStructureCollection_; ///< collection of all available molecular structures 
        Core::Vector startPosition_ = Core::Vector(0.0, 0.0, 0.0); 
        Core::Vector startRotation_ = Core::Vector(0.0, 0.0, 0.0);
    };

}

#endif //IDSIMF_COLLISIONMODEL_MDINTERACTIONS_PRECONSTRUCTED_H
