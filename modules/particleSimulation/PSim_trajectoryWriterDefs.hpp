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
 PSim_trajectoryWriterDefs.hpp

 Convienience definitions for trajectory file writers

 ****************************/

#ifndef Particle_simulation_trajectory_writerdefs
#define Particle_simulation_trajectory_writerdefs

#include <functional>
#include <vector>
#include <string>

//forward declare own classes:
namespace BTree {
    class Tree;
    class Particle;
}


namespace ParticleSimulation {

    //function to transform additional parameters
    typedef std::function
            <std::vector<double>(BTree::Particle* particle)>
            additionalPartParamFctType;

    //a name value pair representing an additional parameter to be exported:
    typedef std::pair<std::string, double> additionalParamPair;

    //a empty additional transform function which returns nothing but an empty vector
    static const additionalPartParamFctType emptyParameterTransformFct =
            [](BTree::Particle *) -> std::vector<double>{ return std::vector<double>(); };
    static const std::vector <additionalParamPair> emptyTimestepAdditionalParameters;
}

#endif //Particle_simulation_trajectory_writerdefs
