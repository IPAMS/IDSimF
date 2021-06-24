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

#include "appUtils_ionDefinitionReading.hpp"
#include "appUtils_simulationConfiguration.hpp"
#include "PSim_ionCloudReader.hpp"
#include "PSim_particleStartZone.hpp"
#include "PSim_boxStartZone.hpp"
#include "PSim_cylinderStartZone.hpp"
#include "PSim_util.hpp"
#include "Core_vector.hpp"

/**
 * Returns if a ion cloud definition file was specified in a simulation configuration
 *
 * @param confRoot a simulation configuration
 */
bool AppUtils::isIonCloudDefinitionPresent(const SimulationConfiguration& simConf) {
    return simConf.isParameter(AppUtils::ION_CLOUD_FILE_KEY);
}


/**
 * Reads the ions defined by an ion cloud definition file specified by a simulation configuration
 * into particles and particle pointer vectors.
 *
 * Note that the ion cloud file location is relative to the given configuration file location (confBasePath)
 *
 * @param particles vector to read the particles to
 * @param particlePtrs vector of particle pointers to read the particle pointers to
 * @param confRoot a simulation configuration
 * @param confBasePath path to the base path of the simulation configuration file with the configuration in confRoot
 */
void AppUtils::readIonDefinitionFromIonCloudFile(
        std::vector<std::unique_ptr<BTree::Particle>>& particles,
        std::vector<BTree::Particle*>& particlePtrs,
        const AppUtils::SimulationConfiguration& simConf) {


    std::string ionCloudFileName = simConf.pathRelativeToConfBasePath(
            simConf.stringParameter("ion_cloud_init_file"));

    ParticleSimulation::IonCloudReader reader = ParticleSimulation::IonCloudReader();
    particles = reader.readIonCloud(ionCloudFileName);
    //prepare a vector of raw pointers
    for (const auto& part : particles){
        particlePtrs.push_back(part.get());
    }
}

/**
 * Gets the random particle start zone defined in an simulation configuration file
 *
 * @param confRoot a simulation configuration
 */
std::unique_ptr<ParticleSimulation::ParticleStartZone> AppUtils::getStartZoneFromIonDefinition(
        const SimulationConfiguration& simConf) {

    std::string ionStartGeom_str = simConf.stringParameter("ion_start_geometry");

    Core::Vector ionsBasePos_m = simConf.vector3dParameter("ion_start_base_position_m");

    std::unique_ptr<ParticleSimulation::ParticleStartZone> particleStartZone;
    if (ionStartGeom_str == "box"){
        Core::Vector ionsBoxSize_m = simConf.vector3dParameter("ion_start_box_size_m");
        particleStartZone = std::make_unique<ParticleSimulation::BoxStartZone>(
                ionsBoxSize_m, ionsBasePos_m);
    }
    else if (ionStartGeom_str == "cylinder"){
        double radius = simConf.doubleParameter("ion_start_radius_m");
        double length = simConf.doubleParameter("ion_start_length_m");
        Core::Vector normal_vector = simConf.vector3dParameter("ion_start_cylinder_normal_vector");
        particleStartZone = std::make_unique<ParticleSimulation::CylinderStartZone>(
                radius, length, normal_vector, ionsBasePos_m);
    }
    else{
        std::stringstream ss;
        ss << "Invalid ion start geometry identifier: " << ionStartGeom_str;
        throw (std::invalid_argument(ss.str()));
    }

    return(particleStartZone);
}


/**
 * Sets the kinetic energy of an ion group according to the kinetic energy definition given
 * in confRoot
 * @param particles particle group to set the kinetic energy for
 * @param confRoot a simulation configuration
 */
void AppUtils::setIonsKineticEnergy(
        std::vector<std::unique_ptr<BTree::Particle>>& particles,
        const SimulationConfiguration& simConf) {

    if (simConf.isParameter("ion_kinetic_energy_eV")){
        double ion_ke = simConf.doubleParameter("ion_kinetic_energy_eV") / Core::JOULE_TO_EV;
        std::vector<double> direction_raw = simConf.doubleVectorParameter("ion_direction_vector");
        Core::Vector ion_dir(direction_raw[0], direction_raw[1], direction_raw[2]);
        Core::Vector ion_dir_normalized(ion_dir * (1.0/ion_dir.magnitude()));
        // iterate through all ion groups
        for (auto &particle: particles){
            double ion_velocity = std::sqrt( 2.0* ion_ke / particle->getMass());
            particle->setVelocity(ion_dir_normalized * ion_velocity);
        }
    }
}

/**
 * Reads a random box or random cylinder ion definition into particles and particle pointer vectors
 *
 * @param particles vector to read the particles to
 * @param particlePtrs vector of particle pointers to read the particle pointers to
 * @param confRoot a simulation configuration
 */
void AppUtils::readRandomIonDefinition(
        std::vector<std::unique_ptr<BTree::Particle>>& particles,
        std::vector<BTree::Particle*>& particlePtrs,
        const SimulationConfiguration& simConf) {

    // ions are not given in an init file, read and init random ion box configuration
    std::vector<unsigned int> nIons = simConf.unsignedIntVectorParameter("n_ions");
    std::vector<double> ionMasses = simConf.doubleVectorParameter("ion_masses");
    std::vector<double> ionCharges = simConf.doubleVectorParameter("ion_charges");
    std::vector<double> ionCollisionDiameters_angstrom = simConf.doubleVectorParameter("ion_collision_gas_diameters_angstrom");

    double ions_tob_range = 0.0;
    if (simConf.isParameter("ion_time_of_birth_range_s")){
        ions_tob_range = simConf.doubleParameter("ion_time_of_birth_range_s");
    }

    std::unique_ptr<ParticleSimulation::ParticleStartZone> particleStartZone = getStartZoneFromIonDefinition(simConf);

    // iterate through all ion groups
    for (std::size_t i = 0; i < nIons.size(); i++) {

        // get ion group parameters
        unsigned int nParticles = nIons.at(i);
        double mass = ionMasses.at(i);
        double charge = ionCharges.at(i);
        double collisionDiameter_m = ionCollisionDiameters_angstrom[i]*1e-10;

        // get actual ions for the group
        std::vector<std::unique_ptr<BTree::Particle>> ions = particleStartZone->getRandomParticlesInStartZone(
                nParticles, charge, ions_tob_range);

        // set particle parameters
        for (std::unique_ptr<BTree::Particle>& ion: ions) {
            ion->setMassAMU(mass);
            ion->setDiameter(collisionDiameter_m);
        }
        setIonsKineticEnergy(ions, simConf);

        // and push particles to vectors containing all particles
        for (std::unique_ptr<BTree::Particle>& ion: ions){
            particlePtrs.push_back(ion.get());
            particles.push_back(std::move(ion));
        }
    }
}

/**
 * Reads a particle definition from a simulation config file into particles and particle pointer vectors
 *
 * @param particles vector to read the particles to
 * @param particlePtrs vector of particle pointers to read the particle pointers to
 * @param confRoot a simulation configuration
 * @param confBasePath path to the base path of the simulation configuration file with the configuration in confRoot
 */
void AppUtils::readIonDefinition(std::vector<std::unique_ptr<BTree::Particle>>& particles,
                                 std::vector<BTree::Particle*>& particlePtrs,
                                 const SimulationConfiguration& simConf) {

    if (AppUtils::isIonCloudDefinitionPresent(simConf)) {
        AppUtils::readIonDefinitionFromIonCloudFile(particles, particlePtrs, simConf);
    }
    else {
        AppUtils::readRandomIonDefinition(particles, particlePtrs, simConf);
    }
}

