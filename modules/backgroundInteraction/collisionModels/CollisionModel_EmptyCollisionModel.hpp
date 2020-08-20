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
 CollisionModel_EmptyCollisionModel.hpp

 Empty collision model: collision model which does nothing

 ****************************/

#ifndef Collision_EmptyModel_hpp
#define Collision_EmptyModel_hpp

#include <stdio.h>
#include "BTree_particle.hpp"
#include "CollisionModel_AbstractCollisionModel.hpp"

namespace CollisionModel{

    class EmptyCollisionModel : public AbstractCollisionModel {

    public:

        void updateModelParameters(BTree::Particle& ion) const override;
        void initializeModelParameters(BTree::Particle& ion) const override;

        void modifyAcceleration(Core::Vector& acceleration,
                                BTree::Particle& ion,
                                double dt) override;
        void modifyVelocity(BTree::Particle& ion,
                            double dt) override;
        void modifyPosition(Core::Vector& position,
                            BTree::Particle& ion,
                            double dt) override;
    };

}


#endif //Collision_EmptyModel_hpp
