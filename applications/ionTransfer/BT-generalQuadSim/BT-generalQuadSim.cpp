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
 BT-generalQuadSim.cpp

 Ion trajectory simulation (including space charge and hard sphere collisions) in an
 quadrupole including gas flow and non ideal geometry given by potential arrays

 ****************************/

#include <iostream>
#include <vector>
#include <ctime>
#include "json.h"
#include "parameterParsing.hpp"
#include "BTree_particle.hpp"
#include "BTree_tree.hpp"
#include "Core_randomGenerators.hpp"
#include "PSim_trajectoryExplorerJSONwriter.hpp"
#include "PSim_interpolatedField.hpp"
#include "PSim_util.hpp"
#include "PSim_verletIntegrator.hpp"
#include "CollisionModel_HardSphere.hpp"

// some constants:
double freq_rf = 1.0e6; //Hz, RF Frequency
double omega_rf = freq_rf * M_PI *2.0;
const double rho_per_pa = 2.504e20; //(particles / m^3) / Pa

int main(int argc, const char * argv[]) {

    // read configuration file ======================================================================
    if (argc <2){
        std::cout << "no conf project name or conf file given"<<std::endl;
        return(0);
    }

    std::string confFileName = argv[1];
    std::string confBasePath = confFileBasePath(confFileName);
    std::cout << confFileName << std::endl;

    Json::Value confRoot = readConfigurationJson(confFileName);
    std::cout<<confRoot<<std::endl;

    std::string projectName = argv[2];
    std::cout << projectName << std::endl;


    // read basic simulation parameters =============================================================
    int timeSteps = intConfParameter("sim_time_steps",confRoot);
    int trajectoryWriteInterval = intConfParameter("trajectory_write_interval", confRoot);
    double dt = doubleConfParameter("dt", confRoot);

    // read interpolated fields ======================
    std::unique_ptr<ParticleSimulation::InterpolatedField> rhoField = readInterpolatedField(confBasePath, "rho_field_file", confRoot);
    std::unique_ptr<ParticleSimulation::InterpolatedField> flowField =readInterpolatedField(confBasePath, "flow_field_file", confRoot);
    std::unique_ptr<ParticleSimulation::InterpolatedField> electricFieldQuadRF = readInterpolatedField(confBasePath, "electric_field_rf_file", confRoot);
    std::unique_ptr<ParticleSimulation::InterpolatedField> electricFieldQuadEntrance = readInterpolatedField(confBasePath, "electric_field_entrance_file", confRoot);

    // read physical and geometrical simulation parameters
    int collisionMode = intConfParameter("collision_mode", confRoot);
    double spaceChargeFactor = doubleConfParameter("space_charge_factor", confRoot);
    double collisionGasMassAmu = doubleConfParameter("collision_gas_mass_amu", confRoot);
    double collisionGasDiameterM = doubleConfParameter("collision_gas_diameter_angstrom", confRoot)*1e-10;
    double backgroundTemperture = doubleConfParameter("background_temperature",confRoot);

    double V_rf = doubleConfParameter("V_rf", confRoot);//600; //volts, RF voltage
    double V_entrance = doubleConfParameter("V_entrance", confRoot);
    double P_factor = doubleConfParameter("P_factor", confRoot);

    double entranceAperture = doubleConfParameter("entrance_aperture_mm", confRoot) / 1000.0;

    double startQLength = doubleConfParameter("start_q_length_mm", confRoot) / 1000.0; //the start position on the q length (x) axis
    double reStartQLength = doubleConfParameter("restart_q_length_mm", confRoot) / 1000.0;

    double maxQLength = doubleConfParameter("max_q_length_mm", confRoot) / 1000.0;
    double qStart = doubleConfParameter("q_start_mm", confRoot) / 1000.0;
    double maxRadius = doubleConfParameter("max_r_mm", confRoot) / 1000.0;

    // read ion configuration ========================
    std::vector<int> nIons = std::vector<int>();
    std::vector<double> ionMasses = std::vector<double>();

    if (confRoot.isMember("n_ions")==true){
        Json::Value n_ions_json = confRoot.get("n_ions",0);
        for (int i=0; i<n_ions_json.size(); i++){
            nIons.push_back(n_ions_json.get(i,0.0).asInt());
        }
    }else{
        throw std::invalid_argument("missing configuration value: n_ions");
    }

    if (confRoot.isMember("ion_masses")==true){
        Json::Value ions_masses_json = confRoot.get("ion_masses",0);
        for (int i=0; i<ions_masses_json.size(); i++){
            ionMasses.push_back(ions_masses_json.get(i,0.0).asDouble());
        }
    }else{
        throw std::invalid_argument("missing configuration value: ions_masses");
    }

    std::vector<std::unique_ptr<BTree::Particle>>particles;
    std::vector<BTree::Particle*>particlePtrs;


    //prepare file writers ==============================================================================
    auto jsonWriter = std::make_unique<ParticleSimulation::TrajectoryExplorerJSONwriter>
            (projectName+ "_trajectories.json");
    jsonWriter->setScales(1000,1e6);

    // prepare random generators:
    Core::RndDistPtr rnd_x = Core::globalRandomGenerator->getUniformDistribution(qStart,reStartQLength);
    Core::RndDistPtr rnd_yz = Core::globalRandomGenerator->getUniformDistribution(-entranceAperture,entranceAperture);

    //init ions:
    for (int i=0; i<nIons.size(); i++){
        int nParticles = nIons[i];
        double mass = ionMasses[i];
        auto ions= ParticleSimulation::util::getRandomIonsInBox(nParticles,1.0,
                                                                   Core::Vector(qStart,-entranceAperture/2.0,-entranceAperture/2.0),
                                                                   Core::Vector(startQLength,entranceAperture,entranceAperture));

        for (int j=0; j<nParticles; j++){
            ions[j]->setMassAMU(mass);
            particlePtrs.push_back(ions[j].get());
            particles.push_back(std::move(ions[j]));
        }
    }


    auto backgroundGasVelocityFunction = [&flowField](Core::Vector& location){
        Core::Vector flowVelo =flowField->getInterpolatedVector(location.x(), location.y(), location.z(), 0);
        return flowVelo;
    };

    auto backgroundGasPressureFunction = [&rhoField,P_factor](Core::Vector& location){
        double rho =rhoField->getInterpolatedScalar(location.x(), location.y(), location.z(), 0);
        double pressure_pa = rho / rho_per_pa * P_factor;

        //std::cout << "rho "<<rho<<" pressure_pa "<<pressure_pa<<std::endl;
        return pressure_pa;
    };


    //init gas collision models:
    CollisionModel::HardSphereModel hsModel = CollisionModel::HardSphereModel(
            backgroundGasPressureFunction,
            backgroundGasVelocityFunction,
            backgroundTemperture,
            collisionGasMassAmu,
            collisionGasDiameterM
            );


    // define functions for the trajectory integration ==================================================
    auto accelerationFunction = [V_rf,V_entrance,spaceChargeFactor,collisionMode,&electricFieldQuadRF,&electricFieldQuadEntrance](BTree::Particle* particle, int particleIndex, BTree::Tree& tree, double time,int timestep){
        //x is the long quad axis
        Core::Vector pos = particle->getLocation();
        double particleCharge =particle->getCharge();

        Core::Vector E =
                (electricFieldQuadRF->getInterpolatedVector(pos.x(), pos.y(), pos.z(), 0) * cos(omega_rf* time)*V_rf) +
                (electricFieldQuadEntrance->getInterpolatedVector(pos.x(), pos.y(), pos.z(), 0)*V_entrance);

        //FIXME: This was intended to test if an ion has left the geometry,
        //with the updated implementations of the interpolated field class, an ion leaving the
        //interpolated field will result in an exception rather than a vector of NaNs.

        if (std::isnan(E.x()) == true){
            if(collisionMode == 1){
                particle->setActive(false);
            } else{
                particle->setInvalid(true);
                return (Core::Vector(0.0,0.0,0.0));
            }
        }

        Core::Vector spaceChargeForce = tree.computeEFieldFromTree(*particle)* spaceChargeFactor;
        Core::Vector result = (E  + spaceChargeForce) * particleCharge / particle->getMass();
        return(result);
    };

    ParticleSimulation::additionalPartParamFctType additionalParameterTransformFct =
            [&backgroundGasPressureFunction](BTree::Particle *particle) -> std::vector<double>{
                double pressure_pa = backgroundGasPressureFunction(particle->getLocation());
                std::vector<double> result = {
                        particle->getVelocity().x(),
                        particle->getVelocity().y(),
                        particle->getVelocity().z(),
                        pressure_pa
                };
                return result;
            };

    auto timestepWriteFunction = [trajectoryWriteInterval, &additionalParameterTransformFct, &jsonWriter](
            std::vector<BTree::Particle *> &particles, BTree::Tree &tree, double time, int timestep, bool lastTimestep){
        if (timestep % trajectoryWriteInterval ==0){

            std::cout<<"ts:"<<timestep<<" time:"<<time<<std::endl;
            char buffer [100];
            int cx;
            cx = snprintf ( buffer, 100, "%06d", timestep );
            jsonWriter->writeTimestep(particles,additionalParameterTransformFct,time,false);
        }
        if (lastTimestep){
            jsonWriter->writeTimestep(particles,additionalParameterTransformFct, time,true);
            jsonWriter->writeSplatTimes(particles);
            jsonWriter->writeIonMasses(particles);
            std::cout<<"finished ts:"<<timestep<<" time:"<<time<<std::endl;
        }
    };

    auto otherActionsFunction= [maxQLength,maxRadius,&rnd_yz,&rnd_x](Core::Vector& newPartPos,BTree::Particle* particle,int particleIndex, BTree::Tree& tree, double time, int timestep){

        double r_pos = std::sqrt( newPartPos.y()*newPartPos.y() + newPartPos.z()*newPartPos.z() );

        if (newPartPos.x() > maxQLength){
            newPartPos.z(0.0);
        }

        if (r_pos > maxRadius || newPartPos.x() > maxQLength || particle->isInvalid() == true){
            double yNew = rnd_yz->rndValue();
            double zNew = rnd_yz->rndValue();

            newPartPos.set(
                    rnd_x->rndValue(),
                    yNew,
                    zNew);
            particle->setInvalid(false);
        }
    };

    // simulate ===============================================================================================
    clock_t begin = std::clock();
    ParticleSimulation::VerletIntegrator verletIntegrator(
            particlePtrs,
            accelerationFunction, timestepWriteFunction, otherActionsFunction,
            hsModel);
    verletIntegrator.run(timeSteps,dt);

    clock_t end = std::clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

    std::cout << particles[0]->getLocation()<<std::endl;
    std::cout << "elapsed secs"<< elapsed_secs<<std::endl;
}