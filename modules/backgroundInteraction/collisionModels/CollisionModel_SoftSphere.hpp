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
 CollisionModel_HardSphere.hpp

 Hard Sphere model for ion / neutral gas particle collisions

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

#ifndef IDSIMF_COLLISIONMODEL_SOFTSPHERE_H
#define IDSIMF_COLLISIONMODEL_SOFTSPHERE_H

#include "CollisionModel_AbstractCollisionModel.hpp"

namespace CollisionModel{

    class SoftSphereModel : public AbstractCollisionModel {

    public:

        void initializeModelParameters(BTree::Particle& ion) const;

        void updateModelParameters(BTree::Particle& ion) const;

        void modifyAcceleration(Core::Vector& acceleration,
                                        BTree::Particle& particle,
                                        double dt);
        void modifyVelocity(BTree::Particle& particle,
                                    double dt);
        void modifyPosition(Core::Vector& position,
                                    BTree::Particle& particle,
                                    double dt);
    };

}

#endif //IDSIMF_COLLISIONMODEL_SOFTSPHERE_H
