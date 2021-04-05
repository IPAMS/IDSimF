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

#include "PSim_trajectoryExplorerJSONwriter.hpp"
#include "BTree_tree.hpp"
#include "BTree_particle.hpp"

/**
 * Constructor: Creates the JSON file writer and inits the JSON file to write to
 * @param jsonFilename a filename of the file to write to
 */
ParticleSimulation::TrajectoryExplorerJSONwriter::TrajectoryExplorerJSONwriter(std::string jsonFilename){
    jsonFile_ = std::make_unique<std::ofstream>();
    jsonFile_->open(jsonFilename);
    initFile();
}

/**
 * Destructor: Closes the JSON file and destroys the file writer
 */
ParticleSimulation::TrajectoryExplorerJSONwriter::~TrajectoryExplorerJSONwriter(){
    closeFile();
}

/**
 * Sets the spatial and temporal scaling factors
 *
 * @param scale spatial scaling factor
 * @param timeScale temporal scaling factor
 */
void ParticleSimulation::TrajectoryExplorerJSONwriter::setScales(double scale, double timeScale){
    this->scale_ = scale;
    this->timeScale_ = timeScale;
}

/**
 * Writes the current timestep to the JSON file without additional parameters.
 *
 * \see ParticleSimulation::TrajectoryExplorerJSONwriter::writeTimestep(Core::Tree&,std::vector<int>,
 *      const ParticleSimulation::partAttribTransformFctType&,
 *      const std::vector<additionalParamPair>,
 *      double,
 *      bool)
 *
 * for details
 */
void ParticleSimulation::TrajectoryExplorerJSONwriter::writeTimestep(std::vector<BTree::Particle*>& particles,
                                                                     double time, bool lastTime){
    writeTimestep(particles, emptyParameterTransformFct, emptyTimestepAdditionalParameters, time,
                  lastTime);
}

/**
 * Writes the current timestep to the JSON file with additional parameters for the particles.
 *
 * \see ParticleSimulation::TrajectoryExplorerJSONwriter::writeTimestep(Core::Tree&,std::vector<int>,
 *      const ParticleSimulation::partAttribTransformFctType&,
 *      const std::vector<additionalParamPair>,
 *      double,
 *      bool)
 *
 * for details
 */
void ParticleSimulation::TrajectoryExplorerJSONwriter::writeTimestep(
        std::vector<BTree::Particle*>& particles,
        const partAttribTransformFctType &particleParameterTransformFct,
        double time,
        bool lastTime){

    writeTimestep(particles,particleParameterTransformFct,emptyTimestepAdditionalParameters,time,lastTime);
}

/**
 * Writes the current state of a Core, which represents the state of a timestep, to the JSON file.
 * The particles which are written to the JSON file are specified / filtered out by their "external keys" which
 * they are associated with in Core.
 *
 * An arbitrary set of additional parameters can be exported to the JSON file for the current time step.
 * An additional parameter is defined by a name and value pair (string / double) pair.
 *
 * The particle positions are always written to the time step, an arbitrary function is
 * used to transform additional particle parameters to output written in the JSON file.
 * This function returns a vector of double values and takes a single particle reference.
 * (If no additional data should be written, use emptyParameterTransformFct which only returns empty vectors)
 *
 * @param tree a tree which is written to the JSON file
 * @param externalParticleKeys external keys of the particles which are written to the JSON file
 * @param particleParameterTransformFct transform function to transform a particle into numeric additional parameters
 * @param timestepAdditionalParameters vector of additional parameter pairs for this time step
 * @param time the time of the current time step
 * @param lastTime must be true if this is the last time step to be written to the JSON file
 */
void ParticleSimulation::TrajectoryExplorerJSONwriter::writeTimestep(
        std::vector<BTree::Particle*>& particles,
        const ParticleSimulation::partAttribTransformFctType &particleParameterTransformFct,
        const std::vector<additionalParamPair> timestepAdditionalParameters,
        double time,
        bool lastTime){

    *jsonFile_<<"{\"time\":"<<time*this->timeScale_;

    for (auto& paramPair: timestepAdditionalParameters){
        *jsonFile_ << ",\"" << paramPair.first <<"\":"<<paramPair.second;
    }

    *jsonFile_<<",\n"<<"\"ions\":\n[";

    for(auto it = particles.begin(); it != particles.end(); ++it){
        BTree::Particle *particle = *it;
        std::vector<double> additionalParameters = particleParameterTransformFct(particle);

        if (! additionalParameters.empty()){
            *jsonFile_ << "[";
            this->writeIonPosition(particle);
            for (auto &val : additionalParameters) {
                *jsonFile_ << "," << val;
            }
            *jsonFile_ << "]";
        }
        else{
            this->writeIonPosition(particle);
        }

        if(std::next(it) != particles.end()) { //not the last element: write comma
            *jsonFile_ << ",\n";
        }
    }

    if (! lastTime) {
        *jsonFile_ << "]},\n";
    } else {
        *jsonFile_ << "]}\n]";
    }
}

/**
 * Writes the particle splat times (times when the particles have terminated) to the JSON file
 *
 * @param tree a tree containing particles which splat times are written to the JSON file
 * @param externalParticleKeys external keys of the particles which splat times are written to the JSON file
 */
void ParticleSimulation::TrajectoryExplorerJSONwriter::writeSplatTimes(std::vector<BTree::Particle*>& particles){

    *jsonFile_ << ",\n\"splatTimes\":[\n";

     for(auto it = particles.begin(); it != particles.end(); ++it){
         if(std::next(it) != particles.end()) {
            this->writeIonSplatTime(*it, false);
        } else {
            this->writeIonSplatTime(*it, true);
        }
    }
    *jsonFile_ << "]";
}

/**
 * Writes the particle masses to the JSON file
 *
 * @param tree a tree containing particles which masses are written to the JSON file
 * @param externalParticleKeys external keys of the particles which masses are written to the JSON file
 */
void ParticleSimulation::TrajectoryExplorerJSONwriter::writeIonMasses(std::vector<BTree::Particle*>& particles){

    *jsonFile_ << ",\n\"ionMasses\":[\n";

     for(auto it = particles.begin(); it != particles.end(); ++it){

        *jsonFile_ << (*it)->getMass() / Core::AMU_TO_KG;
         if(std::next(it) != particles.end()) {
            *jsonFile_ << ",\n";
        } else {
            *jsonFile_ << "\n";
        }
    }
    *jsonFile_ << "]\n";
}

/**
 * Initializes the JSON file
 */
void ParticleSimulation::TrajectoryExplorerJSONwriter::initFile(){
    std::string header = "{\n"
                         "\"steps\":[";
    *jsonFile_ << header;
}

/**
 * Closes the JSON file
 */
void ParticleSimulation::TrajectoryExplorerJSONwriter::closeFile(){
    std::string footer = "}";
    *jsonFile_ << footer;
    jsonFile_->flush();
}

/**
 * Writes the position of a particle to the JSON file
 *
 * @param particle a particle which position should be written to the JSON file
 */
void ParticleSimulation::TrajectoryExplorerJSONwriter::writeIonPosition(BTree::Particle* particle){

    const Core::Vector* ionPos = &particle->getLocation();
    *jsonFile_ << "[" << ionPos->x() * scale_ << "," << ionPos->y() * scale_ << "," << ionPos->z() * scale_ << "]";

}

/**
 * Write the splat time of a particle to the JSON file
 *
 * @param particle a particle which splat time should be written to JSON file
 * @param lastParticle has to be true for the last splat time written to the JSON file
 * (to close the vector of splat times in the JSON file)
 */
void ParticleSimulation::TrajectoryExplorerJSONwriter::writeIonSplatTime(BTree::Particle *particle, bool lastParticle){

    *jsonFile_ << particle->getSplatTime() * this->timeScale_;
    if (lastParticle) {
        *jsonFile_ << "\n";
    } else {
        *jsonFile_ << ",\n";
    }
}