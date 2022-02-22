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

 ------------
 PSim_util.hpp

 Some convenience functions for particle sampling / particle generation

 ****************************/

#ifndef BTree_particleSimulation_
#define BTree_particleSimulation_

#include <vector>
#include <tuple>
#include <memory>

//forward declare own classes:
namespace Core{
    class Vector;
}
namespace BTree{
    class Particle;
    class Tree;
}

namespace ParticleSimulation {
    enum Plane {XY,XZ,YZ};

    namespace util {
        [[nodiscard]] std::vector<std::unique_ptr<BTree::Particle>>
            prepareIonsOnCylinderWalls(unsigned int nIons, double charge, double radius, double length);

        [[nodiscard]] std::vector<std::tuple<double,double,Core::Vector>>
            probeForces(std::vector<BTree::Particle>& ions, Plane plane, int nU, int nV, double minU, double minV,
                    double maxU, double maxV, double slicePos);

        [[nodiscard]] std::vector<Core::Vector> getRandomPositionsInBox(unsigned int nPositions, Core::Vector corner, Core::Vector boxSize);

        //std::vector<std::unique_ptr<BTree::Particle>>
        //getRandomIonsInBox(int numIons, double charge, Core::Vector corner, Core::Vector boxSize, double timeOfBirthRange=0.0);

        [[nodiscard]] std::vector<std::unique_ptr<BTree::Particle>>
            getIonOnLineVector(unsigned int numIons, double charge, double x, double y, double z, double timeOfBirthRange=0.0);
    }
}

#endif
