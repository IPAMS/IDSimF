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

#ifndef Collision_HardSphere_hpp
#define Collision_HardSphere_hpp

#include "BTree_particle.hpp"
#include "Core_constants.hpp"
#include "CollisionModel_AbstractCollisionModel.hpp"
#include "CollisionModel_SpatialFieldFunctions.hpp"
#include "CollisionModel_MathFunctions.hpp"
#include "RS_AbstractReaction.hpp"
#include <stdio.h>
#include <functional>

namespace CollisionModel{

    class HardSphereModel : public AbstractCollisionModel  {

    public:
        constexpr static double DIAMETER_N2 = 3.64e-10;
        constexpr static double DIAMETER_HE = 2.80e-10;

        HardSphereModel(
                double staticPressure,
                double staticTemperature,
                double collisionGasMassAmu,
                double collisionGasDiameterM,
                bool maxwellianApproximation = false);

        HardSphereModel(
                std::function<double(Core::Vector& location)>pressureFunction,
                std::function<Core::Vector(Core::Vector& location)>velocityFunction,
                double StaticTemperature,
                double collisionGasMassAmu,
                double collisionGasDiameterM,
                bool maxwellianApproximation = false);

        HardSphereModel(
                double staticPressure,
                double staticTemperature,
                double collisionGasMassAmu,
                double collisionGasDiameterM,
                std::function<void(RS::CollisionConditions, BTree::Particle&)>afterCollisionFunction,
                bool maxwellianApproximation = false);

        HardSphereModel(
                std::function<double(Core::Vector& location)>pressureFunction,
                std::function<Core::Vector(Core::Vector& location)>velocityFunction,
                double StaticTemperature,
                double collisionGasMassAmu,
                double collisionGasDiameterM,
                std::function<void(RS::CollisionConditions, BTree::Particle&)>afterCollisionFunction,
                bool maxwellianApproximation = false);


        void updateModelParameters(BTree::Particle& ion) const override;
        void initializeModelParameters(BTree::Particle& ion) const override;

        void modifyAcceleration(
                Core::Vector& acceleration,
                BTree::Particle& ion,
                double dt) override;

        void modifyVelocity(
                BTree::Particle& ion,
                double dt) override;

        void modifyPosition(
                Core::Vector& position,
                BTree::Particle& ion,
                double dt) override;

    private:
        bool maxwellianApproximation_ = false;  ///< flag if a pure maxwellian approximation for the gas particles is used

        double temperature_K_ = 298;    ///< static background temperature in K
        double collisionGasMass_Amu_;   ///< mass of the neutral colliding gas particles in amu
        double collisionGasMass_kg_;    ///< mass of the neutral colliding gas particles in kg
        double collisionGasDiameter_m_; ///< effective collision diameter of the neutral collision gas particles in m

        const double PI_SQRT = std::sqrt(M_PI);
        const double PI_2 = 2.0*M_PI;
        const double SQRT3_3 = std::sqrt(3) * 3;

        std::function<double(Core::Vector&)> pressureFunction_ = nullptr; ///< a spatial pressure function
        std::function<Core::Vector(Core::Vector&)> velocityFunction_ = nullptr; ///< a spatial velocity function
        std::function<void(RS::CollisionConditions, BTree::Particle&)> afterCollisionActionFunction_ = nullptr;
        ///< Function with things to do after a collision (e.g. collision based chemical reactions)
    };
}

#endif /* Collision_HardSphere_hpp */
