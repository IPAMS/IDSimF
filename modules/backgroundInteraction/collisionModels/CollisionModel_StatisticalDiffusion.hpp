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
 CollisionModel_StatisticalDiffusion.hpp

 C++ implementation of the statistical diffusion model
 described by Appelhans and Dahl in
 "SIMION ion optics simulations at atmospheric pressure"
 DOI: 10.1016/j.ijms.2005.03.010

 ****************************/

#ifndef Collision_StatisticalDiffusion_hpp
#define Collision_StatisticalDiffusion_hpp

#include <stdio.h>
#include "Core_vector.hpp"
#include "BTree_particle.hpp"
#include "Core_constants.hpp"
#include "CollisionModel_util.hpp"
#include "CollisionModel_AbstractCollisionModel.hpp"
#include "CollisionModel_SpatialFieldFunctions.hpp"
#include "CollisionModel_MathFunctions.hpp"
#include "CollisionModel_CollisionStatistics.hpp"

namespace CollisionModel {
    class StatisticalDiffusionModel : public AbstractCollisionModel {
        public:

            constexpr static int index_ionSTPdamping = 0;
            constexpr static int index_tRatio  = 1;
            constexpr static int index_ptRatio  = 2;

            StatisticalDiffusionModel(
                    double staticPressure,
                    double staticTemperature,
                    double collisionGasMassAmu,
                    double collisionGasDiameterM,
                    CollisionStatistics cs = CollisionStatistics());

            StatisticalDiffusionModel(
                    double staticPressure,
                    double staticTemperature,
                    Core::Vector staticGasVelocity,
                    double collisionGasMassAmu,
                    double collisionGasDiameterM,
                    CollisionStatistics cs = CollisionStatistics());

            StatisticalDiffusionModel(
                    std::function<double(const Core::Vector& location)>pressureFunction,
                    std::function<double(const Core::Vector& location)>temperatureFunction,
                    std::function<Core::Vector(const Core::Vector& location)>velocityFunction,
                    double collisionGasMassAmu,
                    double collisionGasDiameterM,
                    CollisionStatistics cs = CollisionStatistics());

            void setSTPParameters(BTree::Particle& ion) const;
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

        private:
            constexpr double static STP_TEMP = 273.15;        ///< Standard temperature (K)
            constexpr double static STP_PRESSURE = 100000.0;  ///< Standard pressure (Pa)
            constexpr double static M2_PER_NM2 = 1.0e-18;     ///< (m^2/nm^2)
            const double STP_PARTICLE_DENSITY =   Core::N_AVOGADRO / Core::MOL_VOLUME; ///< Particles per volume (n/m3) at STP.


        // Spatial parameter functions for the background gas:
            std::function<double(const Core::Vector&)>pressureFunction_;        ///< Spatial pressure function
            std::function<double(const Core::Vector&)>temperatureFunction_;     ///< Spatial temperature function
            std::function<Core::Vector(const Core::Vector&)>velocityFunction_;  ///< Spatial velocity function
            CollisionStatistics cs_; ///< CollisionStatistics
            // FIXME: this vector should be set to const...
            std::vector<std::vector<double>> icdfs_;

            double collisionGasMass_amu_ = 0.0;    ///< Mass of the collision gas particles in amu
            double collisionGasDiameter_nm_ = 0.0; ///< Effective diamenter of collision gas particles in nanometer

            double randomWalkDistance_(double logParticleMassRatio) const;

    };
}


#endif //Collision_StatisticalDiffusion_hpp
