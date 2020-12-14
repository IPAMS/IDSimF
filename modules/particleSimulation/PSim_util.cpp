/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2020 - Physical and Theoretical Chemistry /
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

#include "PSim_util.hpp"
#include "BTree_tree.hpp"
#include "BTree_particle.hpp"
#include "Core_vector.hpp"
#include "Core_randomGenerators.hpp"
#include <random>
#include <cmath>

typedef std::mt19937 MyRNG;  // the Mersenne Twister with a popular choice of parameters

/**
 * Creates a hollow cylinder with ions / charged particles on the cylinder wall.
 * Cylinder axis is the z axis
 *
 * @param nIons number of particles to create
 * @param charge charge of particles to create
 * @param radius radius of the cylinder (in xy direction)
 * @param length length of the cylinder (in z direction)
 * @return particles on the cylinder wall
 */
std::vector<std::unique_ptr<BTree::Particle>> ParticleSimulation::util::prepareIonsOnCylinderWalls(int nIons,double charge, double radius, double length){
    
    MyRNG rng;
    rng.seed(1.0);
    std::uniform_real_distribution<double> rnd_len(-length/2.0,length/2.0);
    std::uniform_real_distribution<double> rnd_phi(0.0,2.0*M_PI);

    std::vector<std::unique_ptr<BTree::Particle>> result;

    for (int i=0; i<nIons; i++){
        double phi = rnd_phi(rng);
        double z = rnd_len(rng);
        double x = radius*sin(phi);
        double y = radius*cos(phi);
        std::unique_ptr<BTree::Particle> newIon = std::make_unique<BTree::Particle>(Core::Vector(x,y,z),charge);
        newIon->setActive(false);
        result.push_back(std::move(newIon));
    }
    
    return(result);
}

/**
 * Probes electrical force induced by a group of charged particles on a spatial plane
 * (U V W are the spatial coordinates in the sampled plane)
 *
 * @param ions cloud of charged particles which induce electric force by their charge
 * @param plane a plane orientation
 * @param nU number of sampled points in U direction of the plane
 * @param nV number of sampled points in V direction of the plane
 * @param minU lower bound of the plane in U direction
 * @param minV lower bound of the plane in V direction
 * @param maxU upper bound of the plane in U direction
 * @param maxV upper bound of the plane in V direction
 * @param slicePos position of the plane in W direction
 * @return a vector with the U and V coordinates and the local electrical forces
 */
std::vector<std::tuple<double,double,Core::Vector>> ParticleSimulation::util::probeForces(
        std::vector<BTree::Particle>& ions,
        ParticleSimulation::Plane plane,
        int nU,
        int nV,
        double minU,
        double minV,
        double maxU,
        double maxV,
        double slicePos){
    
    double dU = (maxU - minU)/nU;
    double dV = (maxV - minV)/nV;
    double uPos;
    double vPos;

    std::vector<std::tuple<double,double,Core::Vector>> result;

    Core::Vector domainMin = Core::Vector(-100,-100,-100);
    Core::Vector domainMax = Core::Vector( 100, 100, 100);
    BTree::Tree tree = BTree::Tree(domainMin,domainMax);
    int k = 0;
    for ( auto &ion : ions ) {
        tree.insertParticle(ion,k);
        ++k;
    }
    tree.computeChargeDistribution();
    
    Core::Vector force(0,0,0);
    BTree::Particle testParticle;
    Core::Vector particlePosition(0,0,0);
    for (int i=0; i<nU; ++i){
        uPos = minU + (i*dU);
        for(int j=0; j<nV; ++j){
            vPos = minV + (j*dV);
            
            if(plane== XY){
                particlePosition = Core::Vector(uPos,vPos,slicePos);
                
            }else if(plane==XZ){
                particlePosition = Core::Vector(uPos,slicePos,vPos);
                
            }else{
                particlePosition = Core::Vector(slicePos,uPos,vPos);
            }
            
            testParticle = BTree::Particle(particlePosition, 1);
            force= tree.computeEFieldFromTree(testParticle);
            result.emplace_back(std::tuple<double,double,Core::Vector>{uPos,vPos,force});
            //resultFile<<uPos<<" "<<vPos<<" "<<force<<std::endl;
        }
    }
    return(result);
}

/**
 * Produces random positions in a box
 * @param nPositions number of positions
 * @param corner lower corner of the box (xLow,yLow,zLow corner)
 * @param boxSize size of the box in x,y,z direction
 * @return random sampled positions in the speciefied box
 */
std::vector<Core::Vector> ParticleSimulation::util::getRandomPositionsInBox(int nPositions, Core::Vector corner,
                                                                             Core::Vector boxSize) {
    std::vector<Core::Vector> result;
    Core::RndDistPtr rnd_x = Core::globalRandomGenerator->getUniformDistribution(0,boxSize.x());
    Core::RndDistPtr rnd_y = Core::globalRandomGenerator->getUniformDistribution(0,boxSize.y());
    Core::RndDistPtr rnd_z = Core::globalRandomGenerator->getUniformDistribution(0,boxSize.z());

    for (int i=0; i<nPositions; i++) {
        result.push_back(
            Core::Vector(
                    corner.x() + rnd_x->rndValue(),
                    corner.y() + rnd_y->rndValue(),
                    corner.z() + rnd_z->rndValue())
        );
    };

    return result;
}

/**
 * Gets randomly sampled ions in a box
 * @param numIons number of ions
 * @param charge charge of the generated ions
 * @param corner lower corner of the box (xLow,yLow,zLow corner)
 * @param boxSize size of the box in x,y,z direction
 * @param timeOfBirthRange ions are generated with times of birth uniformly distributed in this range
 * @return randomly vector of ions in a box
 */
std::vector<std::unique_ptr<BTree::Particle>> ParticleSimulation::util::getRandomIonsInBox(int numIons, double charge,
                                                                                           Core::Vector corner, Core::Vector boxSize,
                                                                                           double timeOfBirthRange){
    std::vector<Core::Vector> positions = getRandomPositionsInBox(numIons,corner,boxSize);
    std::vector<std::unique_ptr<BTree::Particle>> result;
    Core::RndDistPtr rnd_tob = Core::globalRandomGenerator->getUniformDistribution(0,timeOfBirthRange);
    //Core::Particle* result = new Core::Particle[numIons];
    for (int i=0; i<numIons; i++){
        std::unique_ptr<BTree::Particle> newIon = std::make_unique<BTree::Particle>(positions[i],charge);
        newIon -> setTimeOfBirth(rnd_tob->rndValue());
        result.push_back(std::move(newIon));
    }
    return result;
}

/**
 * Generates a starting vector at position x,y,z
 *
 * @param numIons number of ions
 * @param charge charge of the generated ions
 * @param x x position
 * @param y y position
 * @param z z position
 * @param timeOfBirthRange ions are generated with times of birth uniformly distributed in this range
 * @return vector of ions in a straight line
 */
std::vector<std::unique_ptr<BTree::Particle>>
ParticleSimulation::util::getIonOnLineVector(int numIons, double charge,
                                             double x, double y, double z,
                                             double timeOfBirthRange) {
    Core::RndDistPtr rnd_tob = Core::globalRandomGenerator->getUniformDistribution(0, timeOfBirthRange);

    std::vector<std::unique_ptr<BTree::Particle>> result;
    for (int i = 0; i < numIons; i++) {
        std::unique_ptr<BTree::Particle> newIon = std::make_unique<BTree::Particle>(
                Core::Vector(
                        x,
                        y,
                        z
                ),
                charge);

        newIon->setTimeOfBirth(rnd_tob->rndValue());
        result.push_back(std::move(newIon));
    }
    return result;

}
                                          
