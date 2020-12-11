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
#include "PSim_util.hpp"

/**
 * Returns if a ion cloud definition file was specified in a simulation configuration
 *
 * @param confRoot a simulation configuration
 */
bool AppUtils::isIonCloudDefinitionPresent(Json::Value &confRoot) {
    return confRoot.isMember(AppUtils::ION_CLOUD_FILE_KEY);
}


/**
 * Reads the ions defined by an ion cloud definition file specified by a simulation configuration
 * into particles and particle pointer vectors.
 *
 * Note that the ion cloud file location is relative to the given configuration file location (confFilePathStr)
 *
 * @param particles vector to read the particles to
 * @param particlePtrs vector of particle pointers to read the particle pointers to
 * @param confRoot a simulation configuration
 * @param confBasePath path to the base path of the simulation configuration file with the configuration in confRoot
 */
void AppUtils::readIonDefinitionFromIonCloudFile(
        std::vector<std::unique_ptr<BTree::Particle>>& particles,
        std::vector<BTree::Particle*>& particlePtrs,
        Json::Value& confRoot,
        std::string confBasePath) {

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
 * Reads a random box or random cylinder ion definition into particles and particle pointer vectors
 *
 * @param particles vector to read the particles to
 * @param particlePtrs vector of particle pointers to read the particle pointers to
 * @param confRoot a simulation configuration
 */
void AppUtils::readRandomIonDefinition(
        std::vector<std::unique_ptr<BTree::Particle>>& particles,
        std::vector<BTree::Particle*>& particlePtrs,
        Json::Value& confRoot) {

    // ions are not given in an init file, read and init random ion box configuration
    std::vector<int> nIons = intVectorConfParameter("n_ions", confRoot);
    std::vector<double> ionMasses = doubleVectorConfParameter("ion_masses", confRoot);
    std::vector<double> ionCharges = doubleVectorConfParameter("ion_charges",confRoot);
    std::vector<double> ionCollisionDiameters_angstrom = doubleVectorConfParameter("ion_collision_gas_diameters_angstrom", confRoot);


    double ions_tob_range = 0.0;
    if (confRoot.isMember("ion_time_of_birth_range_s")){
        ions_tob_range = doubleConfParameter("ion_time_of_birth_range_s", confRoot);
    }


    std::string ionStartGeom_str = stringConfParameter("ion_start_geometry",confRoot);
    IonRandomStartGeometry ionStartGeom;


    std::vector<double> basePos = doubleVectorConfParameter("ion_start_base_position_m", confRoot);
    Core::Vector ionsBasePos_m(basePos[0], basePos[1], basePos[2]);

    Core::Vector ionStartBoxSize_m;
    Core::Vector ionStartCornerPosition_m;

    double ionStartCylinder_radius;
    double ionStartCylinder_length;

    if (ionStartGeom_str == "box"){
        ionStartGeom = BOX;
        std::vector<double> boxSize = doubleVectorConfParameter("ion_start_box_size_m", confRoot);
        ionStartBoxSize_m = {boxSize[0], boxSize[1], boxSize[2]};
        ionStartCornerPosition_m = Core::Vector(0.0, 0.0, 0.0) - (ionStartBoxSize_m*0.5);
    }
    else if (ionStartGeom_str == "cylinder"){
        ionStartGeom = CYLINDER;
        ionStartCylinder_radius = doubleConfParameter("ion_start_cylinder_radius_m", confRoot);
        ionStartCylinder_length = doubleConfParameter("ion_start_cylinder_length_m", confRoot);
    }
    else{
        std::stringstream ss;
        ss << "Invalid ion start geometry identifier: " << ionStartGeom_str;
        throw (std::invalid_argument(ss.str()));
    }

    for (int i = 0; i < nIons.size(); i++) {
        int nParticles = nIons[i];
        double mass = ionMasses[i];
        double charge = ionCharges[i];
        double collisionDiameter_m = ionCollisionDiameters_angstrom[i]*1e-10;
        std::vector<std::unique_ptr<BTree::Particle>> ions;

        if (ionStartGeom == BOX) {
            ions = ParticleSimulation::util::getRandomIonsInBox(
                    nParticles, charge,
                    ionStartCornerPosition_m,
                    ionStartBoxSize_m,
                    ions_tob_range);

        }
        else if (ionStartGeom == CYLINDER) {
            // FIXME: use cylinder start zone
            //ions = ParticleSimulation::util::getRandomIonsInCylinderXDirection(
            //        nParticles, charge, ionStartCylinder_radius, ionStartCylinder_length, ions_tob_range);
        }

        for (int j = 0; j < nParticles; j++) {
            ions[j]->setMassAMU(mass);
            ions[j]->setDiameter(collisionDiameter_m);

            //shift all ions according to the base pos:
            Core::Vector buf = ions[j]->getLocation() + ionsBasePos_m;
            ions[j]->setLocation(buf);
            particlePtrs.push_back(ions[j].get());
            particles.push_back(std::move(ions[j]));
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
                                 Json::Value& confRoot, std::string confBasePath) {

    if (AppUtils::isIonCloudDefinitionPresent(confRoot)) {
        AppUtils::readIonDefinitionFromIonCloudFile(particles, particlePtrs, confRoot, confBasePath);
    }
    else {
        AppUtils::readRandomIonDefinition(particles, particlePtrs, confRoot);
    }
}

