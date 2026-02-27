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
 MDIntegrator.cpp

 Description

 ****************************/
#include "CollisionModel_MDIntegrator.hpp"

CollisionModel::MDIntegrator::MDIntegrator(std::string collisionMolecule, bool rotationActive,
                                           std::unique_ptr<AbstractMDForceField> forceField):
collisionMolecule_(collisionMolecule),
rotationActive_(rotationActive),
forceField_(std::move(forceField))
{}

/**
 * Leapfrog method to integrate trajectories of particles involved in a collision.
 * The leapfrog method is of second order and symplectic.
 * @param moleculesPtr collection of molecule pointer
 * @param dt timestep length
 * @param finalTime maximum integration time
 * @param requiredRad radius defining the collision sphere, i.e. the distance that needs to be undercut for
 * a collision to be considered (same radius which is used to estimate the collision probability)
 */
bool CollisionModel::MDIntegrator::leapfrogIntern(std::vector<CollisionModel::Molecule*> moleculesPtr, double dt, double finalTime, double requiredRad){

    bool wasHit = false;

    // distances need to be saved so integration can be stopped if particles leave
    // the domain of interest
    std::vector<double> startDistances;
    size_t moleculesPtr_size = moleculesPtr.size();
    for(size_t i = 0; i < moleculesPtr_size; ++i){
        for(size_t j = i+1; j < moleculesPtr_size; ++j){
            startDistances.push_back((moleculesPtr.at(i)->getComPos() - moleculesPtr.at(j)->getComPos()).magnitude());
        }
    }


    int nSteps = int(round(finalTime/dt));

    std::vector<Core::Vector> forceMolecules(moleculesPtr_size);
    std::vector<Core::Vector> torqueMolecules(moleculesPtr_size);

    forceField_->calculateForceField(moleculesPtr, forceMolecules, torqueMolecules);

    // do the first half step for the velocity, as per leapfrog definition

    size_t i = 0;
    for(auto* molecule : moleculesPtr){
        Core::Vector newComVel =  molecule->getComVel() + forceMolecules.at(i) / molecule->getMass() * dt/2;
        molecule->setComVel(newComVel);
        if(rotationActive_){
            Core::Vector newAnglMom = molecule->getAngMom() + torqueMolecules.at(i) * dt/2;
            molecule->setAngMom(newAnglMom);
            Core::Matrix3 rotMatrix = molecule->getRotationMatrix();
            Core::Matrix3 worldInvInertia = rotMatrix*molecule->getInertiaInvMatrix()*rotMatrix.transpose();
            molecule->setAngVel(worldInvInertia*newAnglMom);
        }
        i++;
    }

    // start the actual leapfrog iteration
    for (int j = 0; j < nSteps; j++){

        // time step for the new position
        i = 0;
        for(auto* molecule : moleculesPtr){
            Core::Vector newComPos =  molecule->getComPos() + molecule->getComVel() * dt;
            molecule->setComPos(newComPos);
            if(rotationActive_){
                Core::Matrix3 rotMatrix = molecule->getRotationMatrix();
                Core::Matrix3 worldInvInertia = rotMatrix*molecule->getInertiaInvMatrix()*rotMatrix.transpose();
                Core::Matrix3 inertiaMatrix = molecule->getInertiaMatrix();
                double I1 = inertiaMatrix(0,0), I2 = inertiaMatrix(1,1), I3=inertiaMatrix(2,2);
                Core::Vector omega = worldInvInertia*molecule->getAngMom();
                if(I1 < CollisionModel::Molecule::MININERTIA){
                    omega.x(0.0);
                }
                if(I2 < CollisionModel::Molecule::MININERTIA){
                    omega.y(0.0);
                }
                if(I3 < CollisionModel::Molecule::MININERTIA){
                    omega.z(0.0);
                }

                Core::Matrix3 newRotMatrix = molecule->getRotationMatrix() + CollisionModel::Molecule::calcRotationMatrixUpdate(rotMatrix, omega) * dt;
                molecule->setRotationMatrix(newRotMatrix);
            }
            i++;
        }
        size_t index = 0;
        for(size_t k = 0; k < moleculesPtr_size; ++k){
            for(size_t l = k+1; l < moleculesPtr_size; ++l){
                if((moleculesPtr.at(l)->getComPos() - moleculesPtr.at(k)->getComPos()).magnitude() > startDistances.at(index++)){
                    return wasHit;
                }
                else if((moleculesPtr.at(l)->getComPos() - moleculesPtr.at(k)->getComPos()).magnitude() <= requiredRad){
                    wasHit=true;
                }
            }
        }

        // recalculate the force
        forceField_->calculateForceField(moleculesPtr, forceMolecules,torqueMolecules);
        i = 0;
        // time step for the new velocity
        for(auto* molecule : moleculesPtr){
            Core::Vector newComVel =  molecule->getComVel() + forceMolecules.at(i) / molecule->getMass() * dt;
            molecule->setComVel(newComVel);
            if(rotationActive_){
                Core::Vector newAnglMom = molecule->getAngMom() + torqueMolecules.at(i) * dt;
                molecule->setAngMom(newAnglMom);
                Core::Matrix3 rotMatrix = molecule->getRotationMatrix();
                Core::Matrix3 worldInvInertia = rotMatrix*molecule->getInertiaInvMatrix()*rotMatrix.transpose();
                molecule->setAngVel(worldInvInertia*newAnglMom);
            }
            i++;
        }
        //Write time step to HDF5 trajectory:
        if(hdf5TWriterConf_.recordingActive == true){
            hdf5TrajectoryWriter_->writeTrajectorySample(
                j*dt, dt, *moleculesPtr.at(0), *moleculesPtr.at(1));
        }
    }
    return false;

}

/**
 * Runge-Kutta 4 method to integrate trajectories of particles involved in a collision.
 * The RK4 is of fourth order.
 * This integrator should NOT be used except for testing purposes as the adaptive step size method is faster.
 * @param moleculesPtr collection of molecule pointer
 * @param dt timestep length
 * @param finalTime maximum integration time
 * @param requiredRad radius defining the collision sphere, i.e. the distance that needs to be undercut for
 * a collision to be considered (same radius which is used to estimate the collision probability)
 */
bool CollisionModel::MDIntegrator::rk4Intern(std::vector<CollisionModel::Molecule*> moleculesPtr, double dt, double finalTime,
                                                                    double requiredRad){


    int nSteps = int(round(finalTime/dt));
    size_t nMolecules = moleculesPtr.size();
    std::vector<Core::Vector> forceMolecules(nMolecules);
    std::vector<Core::Vector> torqueMolecules(nMolecules);


    bool wasHit = false;
    std::vector<double> startDistances;
    for(size_t i = 0; i < nMolecules; ++i){
        for(size_t j = i+1; j < nMolecules; ++j){
            startDistances.push_back((moleculesPtr[i]->getComPos() - moleculesPtr[j]->getComPos()).magnitude());
        }
    }

    size_t i = 0;

    for (int j = 0; j < nSteps; j++){

        std::vector<Core::Vector> velocityMolecules(nMolecules);
        std::vector<Core::Vector> positionMolecules(nMolecules);
        i = 0;
        for(auto* molecule : moleculesPtr){
            velocityMolecules.at(i) = molecule->getComVel();
            positionMolecules.at(i) = molecule->getComPos();
            i++;
        }

        std::vector<Core::Vector> initialPositionMolecules(nMolecules);
        std::vector<Core::Vector> initialVelocityMolecules(nMolecules);
        for(size_t k = 0; k < nMolecules; k++){
            initialPositionMolecules.at(k) = Core::Vector( positionMolecules.at(k).x(), positionMolecules.at(k).y(),positionMolecules.at(k).z() );
            initialVelocityMolecules.at(k) = Core::Vector( velocityMolecules.at(k).x(), velocityMolecules.at(k).y(), velocityMolecules.at(k).z() );
        }

        double length[3] = {1./2, 1./2, 1};
        double mass[2];
        i = 0;
        for(auto* molecule : moleculesPtr){
            mass[i] = molecule->getMass();
            i++;
        }
        forceField_->calculateForceField(moleculesPtr, forceMolecules, torqueMolecules);

        std::array<std::array<Core::Vector, 2>, 4> k;
        std::array<std::array<Core::Vector, 2>, 4> l;


        for(size_t q = 0; q < nMolecules; q++){
            k[0][q] = forceMolecules.at(q) * dt / mass[q];
            l[0][q] = velocityMolecules.at(q) * dt;
        }

        for(size_t n = 1; n < 4; n++){
            i = 0;
            for(auto* molecule : moleculesPtr){
                positionMolecules.at(i) = initialPositionMolecules.at(i) + l[n-1][i]*length[n-1];
                molecule->setComPos(positionMolecules.at(i));
                i++;
            }

            forceField_->calculateForceField(moleculesPtr, forceMolecules, torqueMolecules);

            for(i = 0; i < nMolecules; i++){
                k[n][i] = forceMolecules.at(i) * dt / mass[i];
                l[n][i] = (velocityMolecules.at(i) + k[n-1][i]*length[n-1])*dt;
            }

        }

        i = 0;
        for(auto* molecule : moleculesPtr){
            Core::Vector newComPos = initialPositionMolecules.at(i) + (l[0][i]+ l[1][i]*2 + l[2][i]*2 + l[3][i]) * 1./6;
            molecule->setComPos(newComPos);
            Core::Vector newComVel = initialVelocityMolecules.at(i) + (k[0][i]+ k[1][i]*2 + k[2][i]*2 + k[3][i]) * 1./6;
            molecule->setComVel(newComVel);
            i++;
        }

        size_t index = 0;
        for(size_t k = 0; k < nMolecules; ++k){
            for(size_t l = k+1; l < nMolecules; ++l){
                if((moleculesPtr[l]->getComPos() - moleculesPtr[k]->getComPos()).magnitude() > startDistances[index++]){
                    return wasHit;
                }
                if((moleculesPtr[l]->getComPos() - moleculesPtr[k]->getComPos()).magnitude() <= requiredRad){
                    wasHit=true;
                }
            }
        }
    }

    return false;
}

/**
 * Adaptive step size Runge-Kutta-Fehlberg 45 method to integrate trajectories of particles involved in a collision.
 * This method is of fourth order and uses error control of fifth order on the velocity to adaptively
 * increase or decrease the timestep length reducing the overall computation time.
 * @param moleculesPtr collection of molecule pointer
 * @param dt timestep length
 * @param finalTime maximum integration time
 * @param requiredRad radius defining the collision sphere, i.e. the distance that needs to be undercut for
 * a collision to be considered (same radius which is used to estimate the collision probability)
 * @param tolerance defines the allowed error threshold  to control the timestep lengths
 */
bool CollisionModel::MDIntegrator::rk4InternAdaptiveStep(std::vector<CollisionModel::Molecule*> moleculesPtr, double dt, double finalTime,
                                                                    double requiredRad, double tolerance){

    double integrationTimeSum = 0;
    size_t nMolecules = moleculesPtr.size();
    std::vector<Core::Vector> forceMolecules(nMolecules);
    std::vector<Core::Vector> torqueMolecules(nMolecules);

    size_t i = 0;
    int steps = 0;
    double distance = 0.0;

    // distances need to be saved so integration can be stopped if particles leave
    // the domain of interest
    bool wasHit = false;
    std::vector<double> startDistances;
    for(i = 0; i < nMolecules; ++i){
        for(size_t j = i+1; j < nMolecules; ++j){
            startDistances.push_back((moleculesPtr[i]->getComPos() - moleculesPtr[j]->getComPos()).magnitude());
        }
    }

    std::vector<Core::Vector> velocityMolecules(nMolecules);
    std::vector<Core::Vector> positionMolecules(nMolecules);
    std::vector<Core::Vector> initialPositionMolecules(nMolecules);
    std::vector<Core::Vector> initialVelocityMolecules(nMolecules);

    std::vector<Core::Matrix3> rotMolecules(nMolecules);
    std::vector<Core::Vector> angMomentMolecules(nMolecules);
    std::vector<Core::Vector> initialAngMomentMolecules(nMolecules);
    std::vector<Core::Matrix3> initialRotMolecules(nMolecules);

    std::vector<Core::Matrix3> initialInertInvMolecules(nMolecules);
    std::vector<Core::Matrix3> initialInertMolecules(nMolecules);

    double weight[5][6] = {
                            {1./4, 0, 0, 0, 0, 0},
                            {3./32, 9./32, 0, 0, 0, 0},
                            {1932./2197, -7200./2197, 7296./2197, 0, 0, 0},
                            {439./216, -8, 3680./513, -845./4104, 0, 0},
                            {-8./27, 2, -3544./2565, 1859./4104, -11./40, 0}};
    double mass[2];
    std::array<std::array<Core::Vector, 2>, 6> k;
    std::array<std::array<Core::Vector, 2>, 6> l;
    std::array<Core::Vector, 2> newComVelOrder5;
    std::array<Core::Vector, 2> newComPosOrder4;
    std::array<Core::Vector, 2> newComVelOrder4;

    std::array<std::array<Core::Matrix3, 2>, 6> u;
    std::array<std::array<Core::Vector, 2>, 6> v;
    std::array<Core::Vector, 2> newComAngMomOrder4;
    std::array<Core::Vector, 2> newComAngMomOrder5;
    std::array<Core::Matrix3, 2> newComRotOrder4;
    std::array<Core::Matrix3, 2> newComRotOrder5;

    Core::Matrix3 worldInvInertia;

    std::array<std::array<double,2>, 2> R = {0,0,0,0};
    double globalR, globalDelta;
    double pi = 3.14159;
    Core::RandomSource* rndSource = Core::globalRandomGeneratorPool->getThreadRandomSource();

    while(integrationTimeSum < finalTime){

        i = 0;
        for(auto* molecule : moleculesPtr){
            velocityMolecules[i] = molecule->getComVel();
            positionMolecules[i] = molecule->getComPos();
            initialPositionMolecules[i] = molecule->getComPos();
            initialVelocityMolecules[i] = molecule->getComVel();
            mass[i] = molecule->getMass();

            if(rotationActive_){
                angMomentMolecules[i] = molecule->getAngMom();
                rotMolecules[i] = molecule->getRotationMatrix();
                initialAngMomentMolecules[i] = molecule->getAngMom();
                initialRotMolecules[i] = molecule->getRotationMatrix();
                initialInertInvMolecules[i] = molecule->getInertiaInvMatrix();
                initialInertMolecules[i] = molecule->getInertiaMatrix();
            }
            i++;
        }

        forceField_->calculateForceField(moleculesPtr, forceMolecules, torqueMolecules);

        for(size_t q = 0; q < nMolecules; q++){
            k[0][q] = forceMolecules[q] * dt / mass[q];
            l[0][q] = velocityMolecules[q] * dt;
            if(rotationActive_){

                worldInvInertia = rotMolecules[q]*initialInertInvMolecules[q]*rotMolecules[q].transpose();
                Core::Matrix3 inertiaMatrix = initialInertMolecules[q];
                double I1 = inertiaMatrix(0,0), I2 = inertiaMatrix(1,1), I3=inertiaMatrix(2,2);
                Core::Vector omega = worldInvInertia*angMomentMolecules[q];
                if(I1 < CollisionModel::Molecule::MININERTIA){
                    omega.x(0.0);
                }
                if(I2 < CollisionModel::Molecule::MININERTIA){
                    omega.y(0.0);
                }
                if(I3 < CollisionModel::Molecule::MININERTIA){
                    omega.z(0.0);
                }
                //v[0][q] = Core::Vector(0.0, 0.0, 0.0);
                v[0][q] = torqueMolecules[q] * dt;
                u[0][q] = CollisionModel::Molecule::calcRotationMatrixUpdate(rotMolecules[q], omega) *dt;
            }

        }

        for(size_t n = 1; n < 6; n++){
            for(i = 0; i < nMolecules; i++){
                positionMolecules[i] = initialPositionMolecules[i];
                if(rotationActive_){
                    rotMolecules[i] = initialRotMolecules[i];
                }

            }
            for(size_t m = 0; m < 6; m++){
                i = 0;
                for(auto* molecule : moleculesPtr){
                    positionMolecules[i] += l[m][i]*weight[n-1][m];
                    molecule->setComPos(positionMolecules[i]);
                    if(rotationActive_){
                        rotMolecules[i] = rotMolecules[i] + u[m][i]*weight[n-1][m];
                        molecule->setRotationMatrix(rotMolecules[i]);
                    }
                    i++;
                }
            }

            forceField_->calculateForceField(moleculesPtr, forceMolecules, torqueMolecules);

            for(i = 0; i < nMolecules; i++){
                k[n][i] = forceMolecules[i] * dt / mass[i];
                l[n][i] = velocityMolecules[i];

                if(rotationActive_){
                    worldInvInertia = rotMolecules[i]*initialInertInvMolecules[i]*rotMolecules[i].transpose();
                    Core::Matrix3 inertiaMatrix = initialInertMolecules[i];
                    double I1 = inertiaMatrix(0,0), I2 = inertiaMatrix(1,1), I3=inertiaMatrix(2,2);
                    Core::Vector omega = worldInvInertia*angMomentMolecules[i];
                    if(I1 < CollisionModel::Molecule::MININERTIA){
                        omega.x(0.0);
                    }
                    if(I2 < CollisionModel::Molecule::MININERTIA){
                        omega.y(0.0);
                    }
                    if(I3 < CollisionModel::Molecule::MININERTIA){
                        omega.z(0.0);
                    }
                    //v[n][i] = Core::Vector(0.0, 0.0, 0.0);
                    v[n][i] = torqueMolecules[i] * dt;
                    u[n][i] = CollisionModel::Molecule::calcRotationMatrixUpdate(rotMolecules[i], omega);
                }

                for(size_t m = 0; m < 6; m++){
                    l[n][i] += k[m][i]*weight[n-1][m];
                    if(rotationActive_) u[n][i] = u[n][i] +
                            CollisionModel::Molecule::calcRotationMatrixUpdate(rotMolecules[i],worldInvInertia*v[m][i]*weight[n-1][m]);
                }
                l[n][i] = l[n][i]*dt;
                u[n][i] = u[n][i]*dt;
            }

        }

        for(size_t b = 0; b < nMolecules; ++b){
            for(size_t z = b+1; z < nMolecules; ++z){
                distance = (moleculesPtr[z]->getComPos() - moleculesPtr[b]->getComPos()).magnitude();
            }
        }


        for(i = 0; i < 2; i++){
            newComVelOrder5[i] = initialVelocityMolecules[i] + (k[0][i] * 16./135 + k[2][i] * 6656./12825 + k[3][i] * 28561./56430 + k[4][i] * (-9./50) + k[5][i] * 2./55);
            newComPosOrder4[i] = initialPositionMolecules[i] + (l[0][i] * 25./216 + l[2][i] * 1408./2565 + l[3][i] * 2197./4104 + l[4][i] * (-1./5));
            newComVelOrder4[i] = initialVelocityMolecules[i] + (k[0][i] * 25./216 + k[2][i] * 1408./2565 + k[3][i] * 2197./4104 + k[4][i] * (-1./5));
            if(rotationActive_){
                newComAngMomOrder5[i] = initialAngMomentMolecules[i] + (v[0][i] * (16./135) + v[2][i] * (6656./12825) + v[3][i] * (28561./56430) + v[4][i] * (-9./50) + v[5][i] * (2./55));
                newComAngMomOrder4[i] = initialAngMomentMolecules[i] + (v[0][i] * 25./216 + v[2][i] * 1408./2565 + v[3][i] * 2197./4104 + v[4][i] * (-1./5));
                newComRotOrder4[i] = initialRotMolecules[i] + (u[0][i] * (25./216) + u[2][i] * (1408./2565) + u[3][i] * (2197./4104) + u[4][i] * (-1./5));
                newComRotOrder5[i] = initialRotMolecules[i] + (u[0][i] * (16./135) + u[2][i] * (6656./12825) + u[3][i] * (28561./56430) + u[4][i] * (-9./50) + u[5][i] * (2./55));

            }

        }


        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wfloat-equal"
        for(size_t p = 0; p < 2; p++){


            if(fabs(newComVelOrder5[p].magnitude()) != 0)
                R[p][0] = fabs(newComVelOrder4[p].magnitude()-newComVelOrder5[p].magnitude())/fabs(newComVelOrder5[p].magnitude());
            else
                R[p][0] = 1e-15;

            if(rotationActive_){
                if(fabs(Core::norm1(newComRotOrder5[p])) != 0) {
                    R[p][1] = fabs(Core::norm1(newComRotOrder4[p])-Core::norm1(newComRotOrder5[p]))/fabs(Core::norm1(newComRotOrder5[p]));
                }
                else{
                    R[p][1] = 1e-15;
                }
            }
        }
        if(rotationActive_){
            globalR = std::max({R[1][1], std::max({R[0][1],  std::max({R[0][0],R[1][0]}) }) });
        }else{
            globalR = std::max({R[0][0],R[1][0]});
        }

        if (globalR == 0){
            globalR = 1e-15;
        }
        #pragma GCC diagnostic pop

        globalDelta = 0.84 * std::pow((tolerance/globalR), 1./4);
        integrationTimeSum += dt;
        i = 0;
        for(auto* molecule : moleculesPtr){
            molecule->setComPos(newComPosOrder4[i]);
            molecule->setComVel(newComVelOrder4[i]);
            if(rotationActive_){
                molecule->setAngMom(newComAngMomOrder4[i]);
                molecule->setRotationMatrix(newComRotOrder4[i]);
                Core::Matrix3 worldInvI = newComRotOrder4[i]*initialInertInvMolecules[i]*newComRotOrder4[i].transpose();
                molecule->setAngVel(worldInvI*newComAngMomOrder4[i]);
            }
            i++;

        }

        steps++;
        dt = dt * globalDelta;

        //Write time step to HDF5 trajectory:
        if(hdf5TWriterConf_.recordingActive == true){
            hdf5TrajectoryWriter_->writeTrajectorySample(
                integrationTimeSum, dt, *moleculesPtr.at(0), *moleculesPtr.at(1));
        }

        //write time step to legacy trajectory:
        if (legacyTWriterConf_.recordingActive == true) {
            legacyTrajectoryWriter_->writeTrajectorySample(
                integrationTimeSum, dt, moleculesPtr[0]->getComPos(), moleculesPtr[0]->getComVel(),
                moleculesPtr[1]->getComPos(),forceMolecules, distance);
        }

        size_t index = 0;
        for(size_t b = 0; b < nMolecules; ++b){
            for(size_t z = b+1; z < nMolecules; ++z){
                if((moleculesPtr[z]->getComPos() - moleculesPtr[b]->getComPos()).magnitude() > startDistances[index++]){
                    if (legacyTWriterConf_.recordingActive == true) {
                        legacyTrajectoryWriter_->writeTrajectoryDelimiter();
                    }
                    return wasHit;
                }
                if((moleculesPtr[z]->getComPos() - moleculesPtr[b]->getComPos()).magnitude() <= requiredRad){
                    wasHit=true;
                }
            }
        }
    }
    return false;
}


double CollisionModel::MDIntegrator::initRotation(CollisionModel::Molecule& mole, double temperature_K){
    Core::Matrix3 inertiaMatrix = mole.getInertiaMatrix();
    double I1 = inertiaMatrix(0,0), I2 = inertiaMatrix(1,1), I3=inertiaMatrix(2,2);
    double energyRotMolecule = 0;
    if(I1 > 0) {
        energyRotMolecule += 0.5 * Core::K_BOLTZMANN * temperature_K;
    }
    else I1 = 0;
    if(I2 > 0) {
        energyRotMolecule += 0.5 * Core::K_BOLTZMANN * temperature_K;
    }
    else I2 = 0;
    if(I3 > 0) {
        energyRotMolecule += 0.5 * Core::K_BOLTZMANN * temperature_K;
    }
    else I3 = 0;
    double w1 = sqrt(2*energyRotMolecule/(I1+I2+I3));
    Core::Vector anglVelo = {I1 > 0 ? w1 : 0, I2 > 0 ? w1 : 0, I3 > 0 ? w1 : 0};
    mole.setAngVel(anglVelo);
    mole.setAngMom(inertiaMatrix*anglVelo);
    return energyRotMolecule;
}

double CollisionModel::MDIntegrator::calcRotEnergy(Core::Vector omega, Core::Matrix3 I){
    return 0.5*(I(0,0)*omega.x()*omega.x() + I(1,1)*omega.y()*omega.y() + I(2,2)*omega.z()*omega.z());
}


/**
 * Activates trajectory writing and sets trajectory writer configuration
 * @param trajectoryFileName Trajectory Output filename
 * @param trajectoryDistance Distance between ion and background gas in m after which trajectory gets recorded
 * @param startTimeStep First time step which should be written to the trajectory
 * @param minimalSampleInterval Minimal interval between trajectory samples
 */
void CollisionModel::MDIntegrator::setLegacyTrajectoryWriter(const std::string& trajectoryFileName,
                                                              double trajectoryDistance,
                                                              double minimalSampleInterval,
                                                              unsigned int startTimeStep){

    legacyTrajectoryWriter_= std::make_unique<CollisionModel::MDTrajectoryWriter>(trajectoryFileName, minimalSampleInterval);
    legacyTWriterConf_.recordTrajectoryStartTimeStep = startTimeStep;
    legacyTWriterConf_.trajectoryDistance = trajectoryDistance;
    legacyTWriterConf_.modelRecordsTrajectory = true;
}

void CollisionModel::MDIntegrator::setHDF5TrajectoryWriter(const std::string& trajectoryFileName, double trajectoryDistance, unsigned int startTimeStep) {

    hdf5TrajectoryWriter_ = std::make_unique<CollisionModel::HDF5MDTrajectoryWriter>(trajectoryFileName);

    hdf5TWriterConf_.recordTrajectoryStartTimeStep = startTimeStep;
    hdf5TWriterConf_.trajectoryDistance = trajectoryDistance;
    hdf5TWriterConf_.modelRecordsTrajectory = true;
}
