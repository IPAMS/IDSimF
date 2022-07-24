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
 CollisionModel_MultiCollisionModel.hpp

 Aggregate collision model, to combine multiple collision models in one trajectory simulation

 ****************************/

#ifndef Collision_MultiCollisionModel_hpp
#define Collision_MultiCollisionModel_hpp

#include "CollisionModel_AbstractCollisionModel.hpp"

#include <functional>
#include <vector>
#include <memory>

namespace CollisionModel{
    class MultiCollisionModel : public AbstractCollisionModel{

    public:
        explicit MultiCollisionModel(std::vector<std::unique_ptr<AbstractCollisionModel>> models);

        void updateModelParticleParameters(Core::Particle& ion) const override;
        void initializeModelParticleParameters(Core::Particle& ion) const override;
        void updateModelTimestepParameters(int timestep, double time) const override;

        void modifyAcceleration(Core::Vector& acceleration,
                                Core::Particle& ion,
                                double dt) override;
        void modifyVelocity(Core::Particle& ion,
                            double dt) override;
        void modifyPosition(Core::Vector& position,
                            Core::Particle& ion,
                            double dt) override;


    private:
        std::vector<std::unique_ptr<AbstractCollisionModel>> models_;

    };
}


#endif //Collision_MultiCollisionModel_hpp
