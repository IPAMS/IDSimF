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
 ****************************/

#include "CollisionModel_MDInteractionsTrajectorySampler.hpp"
#include "Core_math.hpp"

CollisionModel::MDInteractionsTrajectorySampler::MDInteractionsTrajectorySampler(
            std::unique_ptr<CollisionModel::AbstractMDForceField> forceField,
            bool rotationActive,
            std::unordered_map<std::string, std::shared_ptr<CollisionModel::MolecularStructure>> molecularStructureCollection,
            AppUtils::logger_ptr logger
            ) :
        MDIntegrator("", rotationActive, std::move(forceField)),
        molecularStructureCollection_(std::move(molecularStructureCollection)),
        logger_(logger)
        {}


void CollisionModel::MDInteractionsTrajectorySampler::calculateTrajectory(
        Core::Particle& particle,
        std::string collisionMolecule,
        ParticleInitialConditions moleInitCond, ParticleInitialConditions bgMoleInitCond,
        double integrationTime, double subTimeStep, int maximumSteps, bool ionIsFrozen, MDIntegratorType integratorType) {

    int iterations = 0;

    // Construct the molecule of interest (in most cases the simulated molecular ion) and its atoms
    Molecule mole(moleInitCond.position, moleInitCond.velocity, particle.getMolecularStructure());
    initializeRotation_(mole, moleInitCond.rotationAngles, moleInitCond.angularVelocity);

    // Construct the background gas particle
    Molecule bgMole(bgMoleInitCond.position, bgMoleInitCond.velocity,
                                        molecularStructureCollection_.at(collisionMolecule));
    collisionMolecule_ = collisionMolecule;
    initializeRotation_(bgMole, bgMoleInitCond.rotationAngles, bgMoleInitCond.angularVelocity);

    std::vector<CollisionModel::Molecule*> moleculesPtr = {&mole, &bgMole};

    // possible check for energy conservation
    std::vector<Core::Vector> startVelocity;
    double startEnergy = 0;
    for(auto* molecule : moleculesPtr){
        startVelocity.push_back(molecule->getComVel());
        startEnergy += 0.5 * molecule->getMass() * molecule->getComVel().magnitudeSquared();
    }

    // Call the sub-integrator
    double finalTime = integrationTime; // final integration time in seconds
    double timeStep = subTimeStep; // step size in seconds


    //Start MD integration
    //activate trajectory recording for the whole trajectory:
    if (legacyTWriterConf_.modelRecordsTrajectory) legacyTWriterConf_.recordingActive=true;
    if (hdf5TWriterConf_.modelRecordsTrajectory) {
        hdf5TWriterConf_.recordingActive=true;
        hdf5TrajectoryWriter_->initNewTrajectory(mole.getAtomCount(), bgMole.getAtomCount());
    }

    //trajectorySuccess = rk4Intern(moleculesPtr, timeStep, finalTime, collisionRadius);
    //trajectorySuccess = leapfrogIntern(moleculesPtr, timeStep, finalTime, collisionRadius);
    bool trajectorySuccess;
    // if (integratorType == RK4_ADAPTIVE) {
    //     trajectorySuccess = rk4InternAdaptiveStep(moleculesPtr, timeStep, finalTime, maximumSteps, ionIsFrozen);
    // }
    //std::cout << moleculesPtr[1]->getComPos() << std::endl;
    trajectorySuccess = rk4InternAdaptiveStep(moleculesPtr, timeStep, finalTime, maximumSteps, 100, 1e-8, ionIsFrozen);

    double endEnergy = 0;
    for(auto* molecule : moleculesPtr){
        endEnergy += 0.5 * molecule->getMass() * molecule->getComVel().magnitudeSquared();
    }
    particle.setVelocity(mole.getComVel() + particle.getVelocity());
    ++iterations;
}

void CollisionModel::MDInteractionsTrajectorySampler::initializeRotation_(
        CollisionModel::Molecule& mole, Core::Vector rotationAngles, Core::Vector angularVelocity) {

    Core::Matrix3 initialRotationMatrix = mole.calcRotationMatrix(
    rotationAngles.x(), rotationAngles.y(), rotationAngles.z()
    );

    mole.setRotationMatrix(initialRotationMatrix);
    mole.rotateMoleculeRotationMatrix();

    Core::Matrix3 inertiaMatrix = mole.getInertiaMatrix();
    mole.setAngVel(angularVelocity);
    mole.setAngMom(inertiaMatrix*angularVelocity);
}
