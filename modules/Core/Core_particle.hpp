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

 BTree_particle.hpp

 Simulated charged particle

 ****************************/

#ifndef IDSIMF_CORE_PARTICLE_HPP
#define IDSIMF_CORE_PARTICLE_HPP

#include "Core_constants.hpp"
#include "Core_vector.hpp"
#include <unordered_map>
#include <array>
#include <memory>

//forward declare molecular collision model classes:
namespace CollisionModel{
    class MolecularStructure;
}

namespace Core {

    /**
     * Defines a simulated particle, which can be managed by a geometrical Core
     */
    class Particle {

    public:

        //Constructors:
        Particle() = default;
        // Particle(const Particle& part);
        // Particle(Particle&& part);
        virtual ~Particle()= default;
        Particle(const Core::Vector &location, double chargeElemCharges);
        Particle(const Core::Vector &location, const Core::Vector &velocity, double chargeElemCharges, double massAMU);
        Particle(const Core::Vector &location, const Core::Vector &velocity, double chargeElemCharges, double massAMU, double timeOfBirth);
        Particle(const Core::Vector &location, const Core::Vector &velocity, double chargeElemCharges, double massAMU, double collisionDiameterM, double timeOfBirth);
        Particle(const Core::Vector &location, const Core::Vector &velocity, double chargeElemCharges, 
                    double massAMU, double collisionDiameterM, double timeOfBirth, std::shared_ptr<CollisionModel::MolecularStructure> moleculeStructure);

        // Core::Particle& operator=(const Core::Particle& part);
        // Core::Particle& operator=(Core::Particle&& part);


        //setters and getters for the members of the particle:
        void setLocation(Core::Vector location);
        [[nodiscard]] Core::Vector& getLocation();
        void setVelocity(Core::Vector velocity);
        [[nodiscard]] Core::Vector& getVelocity();
        void setAcceleration(Core::Vector velocity);
        [[nodiscard]] Core::Vector& getAcceleration();

        void setIndex(size_t index);
        [[nodiscard]] std::size_t getIndex() const;

        void setChargeElementary(double chargeElemCharges);
        [[nodiscard]] double getCharge() const;
        void setActive(bool active);
        [[nodiscard]] bool isActive() const;
        void setInvalid(bool invalid);
        [[nodiscard]] bool isInvalid() const;

        [[nodiscard]] double getFloatAttribute(const std::string& key) const;
        void setFloatAttribute(const std::string& key, double value);
        [[nodiscard]] int getIntegerAttribute(const std::string& key) const;
        void setIntegerAttribute(const std::string& key, int value);
        std::array<double, 3>& getAuxCollisionParams();

        void setMobility(double mobility);
        [[nodiscard]] double getMobility() const;
        void setLowFieldMobility(double lowFieldMobility);
        [[nodiscard]] double getLowFieldMobility() const;
        void setMeanFreePathSTP(double meanFreePathSTP);
        [[nodiscard]] double getMeanFreePathSTP() const;
        void setMeanThermalVelocitySTP(double meanVelocitySTP);
        [[nodiscard]] double getMeanThermalVelocitySTP() const;
        void setMassAMU(double massAMU);
        [[nodiscard]] double getMass() const;
        void setDiameter(double diameter);
        [[nodiscard]] double getDiameter() const;

        void setTimeOfBirth(double timeOfBirth);
        [[nodiscard]] double getTimeOfBirth() const;
        void setSplatTime(double splatTime);
        [[nodiscard]] double getSplatTime() const;

        void setMolecularStructure(std::shared_ptr<CollisionModel::MolecularStructure> molecularStructurePtr);
        [[nodiscard]] std::shared_ptr<CollisionModel::MolecularStructure> getMolecularStructure();

    private:
        //all internal variable are in SI units
        std::size_t index_ = 0; ///< an external index (mostly the index of the particle in the reaction model)
        Core::Vector location_= Core::Vector(0.0, 0.0, 0.0); ///< The spatial location of the particle (m)
        Core::Vector velocity_ = Core::Vector(0.0, 0.0, 0.0); ///< The velocity of the particle (m/s)
        Core::Vector acceleration_ = Core::Vector(0.0, 0.0, 0.0); ///< The acceleration of the particle (m/s^2)
        double charge_ = 0.0;           ///< The charge of the particle (C)
        double lowFieldMobility_ = 0.0; ///< The low field (with vanishing electric field) electrical mobility of the charged particle (m^2 / (V*s))
        double mobility_ = 0.0;         ///< The actual electrical mobility of the charged particle (m^2 / (V*s))
        double mass_ = 0.0;             ///< The mass of the particle (kg)
        double diameter_ = 0.0;         ///< (collision) diameter of the particle (m^2)
        double STP_meanFreePath_ = 0.0; ///< mean free path at standard temperature and pressure (m)
        double STP_meanThermalVelocity_ = 0.0; ///< mean thermal velocity at standard temperature and pressure (m/s)

        bool active_ = true;       ///< Flag if the particle is active in the simulation, false if the particle is already terminated
        bool invalid_ = false;     ///< Flag if the particle has an invalid state in the simulation
        double timeOfBirth_ = 0.0; ///< Time when the particle was created
        double splatTime_= 0.0;    ///< Time when the particle was terminated ("splatted")

        std::unordered_map<std::string, double> attributesFloat_; ///< an arbitrary set of additional floating point attributes, accessible by a name
        std::unordered_map<std::string, int> attributesInteger_; ///< an arbitrary set of additional integer attributes, accessible by a name
        std::array<double, 3> auxCollisionParams_ {0.0,0.0,0.0}; ///< quickly accessible parameters for collision models
        std::shared_ptr<CollisionModel::MolecularStructure> molstrPtr;
    };

    typedef std::unique_ptr<Core::Particle> uniquePartPtr;
}


#endif /* IDSIMF_CORE_PARTICLE_HPP */