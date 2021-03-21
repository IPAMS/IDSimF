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
#include "appUtils_parameterParsing.hpp"
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
bool AppUtils::isIonCloudDefinitionPresent(const Json::Value& confRoot) {
    return confRoot.isMember(AppUtils::ION_CLOUD_FILE_KEY);
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
        const Json::Value& confRoot,
        const std::string& confBasePath) {

    std::string ionCloudFileName = pathRelativeToConfBasePath(
            confBasePath,
            confRoot.get("ion_cloud_init_file", 0).asString());

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
std::unique_ptr<ParticleSimulation::ParticleStartZone> AppUtils::getStartZoneFromIonDefinition(const Json::Value &confRoot) {

    std::string ionStartGeom_str = stringConfParameter("ion_start_geometry",confRoot);

    Core::Vector ionsBasePos_m =  vector3dConfParameter("ion_start_base_position_m", confRoot);

    std::unique_ptr<ParticleSimulation::ParticleStartZone> particleStartZone;
    if (ionStartGeom_str == "box"){
        Core::Vector ionsBoxSize_m = vector3dConfParameter("ion_start_box_size_m", confRoot);
        particleStartZone = std::make_unique<ParticleSimulation::BoxStartZone>(
                ionsBoxSize_m, ionsBasePos_m);
    }
    else if (ionStartGeom_str == "cylinder"){
        double radius = doubleConfParameter("ion_start_radius_m", confRoot);
        double length = doubleConfParameter("ion_start_length_m", confRoot);
        Core::Vector normal_vector = vector3dConfParameter("ion_start_cylinder_normal_vector", confRoot);
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
        const Json::Value& confRoot) {

    if (confRoot.isMember("ion_kinetic_energy_eV")){
        double ion_ke = doubleConfParameter("ion_kinetic_energy_eV", confRoot) / Core::JOULE_TO_EV;
        std::vector<double> direction_raw = doubleVectorConfParameter("ion_direction_vector", confRoot);
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
        const Json::Value& confRoot) {

    // ions are not given in an init file, read and init random ion box configuration
    std::vector<int> nIons = intVectorConfParameter("n_ions", confRoot);
    std::vector<double> ionMasses = doubleVectorConfParameter("ion_masses", confRoot);
    std::vector<double> ionCharges = doubleVectorConfParameter("ion_charges",confRoot);
    std::vector<double> ionCollisionDiameters_angstrom = doubleVectorConfParameter("ion_collision_gas_diameters_angstrom", confRoot);

    double ions_tob_range = 0.0;
    if (confRoot.isMember("ion_time_of_birth_range_s")){
        ions_tob_range = doubleConfParameter("ion_time_of_birth_range_s", confRoot);
    }

    std::unique_ptr<ParticleSimulation::ParticleStartZone> particleStartZone = getStartZoneFromIonDefinition(confRoot);

    // iterate through all ion groups
    for (int i = 0; i < nIons.size(); i++) {

        // get ion group parameters
        int nParticles = nIons[i];
        double mass = ionMasses[i];
        double charge = ionCharges[i];
        double collisionDiameter_m = ionCollisionDiameters_angstrom[i]*1e-10;

        // get actual ions for the group
        std::vector<std::unique_ptr<BTree::Particle>> ions = particleStartZone->getRandomParticlesInStartZone(
                nParticles, charge, ions_tob_range);

        // set particle parameters
        for (std::unique_ptr<BTree::Particle>& ion: ions) {
            ion->setMassAMU(mass);
            ion->setDiameter(collisionDiameter_m);
        }
        setIonsKineticEnergy(ions, confRoot);

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
                                 const Json::Value& confRoot,
                                 const std::string& confBasePath) {

    if (AppUtils::isIonCloudDefinitionPresent(confRoot)) {
        AppUtils::readIonDefinitionFromIonCloudFile(particles, particlePtrs, confRoot, confBasePath);
    }
    else {
        AppUtils::readRandomIonDefinition(particles, particlePtrs, confRoot);
    }
}

