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
 idealizedQuadSim.cpp

 Ion trajectory simulation (including space charge and hard sphere collisions) in an idealized
 quadrupole

 ****************************/


#include "Core_particle.hpp"
#include "BTree_tree.hpp"
#include "FileIO_simpleVTKwriter.hpp"
#include "FileIO_trajectoryExplorerJSONwriter.hpp"
#include "PSim_util.hpp"
#include "Integration_verletIntegrator.hpp"
#include "CollisionModel_HardSphere.hpp"
#include "parameterParsing_macro.hpp"
#include "appUtils_stopwatch.hpp"
#include "json.h"
#include <iostream>
#include <vector>
#include <ctime>


int timeSteps = 0;
int trajectoryWriteInterval = 0;
std::unique_ptr<FileIO::TrajectoryExplorerJSONwriter> jsonWriter;

// some constants:

double rho = 8.38/1000.0; //m, inner inscribed diameter
int N_rods = 2; //2N_rods = number of multipole rods
double quadFactor = 0.0;
double quadExp = 2*N_rods-3;

int main(int argc, const char * argv[]) {

    // read configuration file ======================================================================
    if (argc <2){
        std::cout << "no conf project name or conf file given"<<std::endl;
        return(0);
    }

    std::string projectName = argv[1];
    std::cout << projectName<<std::endl;

    std::string confFileName = argv[2];
    std::cout << confFileName<<std::endl;

    std::ifstream confFile;
    confFile.open(confFileName);
    Json::Value confRoot;
    confFile>>confRoot;
    confFile.close();
    std::cout<<confRoot<<std::endl;

    // read basic simulation parameters =============================================================

    double V_rf = 0.0;//600; //volts, RF voltage
    double freq_rf = 1.228e6; //Hz, RF Frequency
    double omega_rf = freq_rf * M_PI *2.0;
    int simType = 0;
    double dt = 0.0;
    std::vector<int> nIons = std::vector<int>();
    std::vector<double> ionMasses = std::vector<double>();

    jsonWriter = std::unique_ptr<FileIO::TrajectoryExplorerJSONwriter>(
            new FileIO::TrajectoryExplorerJSONwriter(projectName+ "_trajectories.json")
    );
    jsonWriter->setScales(1000,1e6);

    intConfParameter("sim_time_steps", timeSteps);
    intConfParameter("trajectory_write_interval", trajectoryWriteInterval);
    intConfParameter("sim_type",simType);
    doubleConfParameter("dt", dt);


    if (confRoot.isMember("n_ions")==true){
        Json::Value n_ions_json = confRoot.get("n_ions",0);
        for (int i=0; i<n_ions_json.size(); i++){
            nIons.push_back(n_ions_json.get(i,0.0).asInt());
        }
    }else{
        std::cout << "missing configuration value: n_ions"<<std::endl;
        return(0);
    }

    if (confRoot.isMember("ion_masses")==true){
        Json::Value ions_masses_json = confRoot.get("ion_masses",0);
        for (int i=0; i<ions_masses_json.size(); i++){
            ionMasses.push_back(ions_masses_json.get(i,0.0).asDouble());
        }
    }else{
        std::cout << "missing configuration value: ions_masses"<<std::endl;
        return(0);
    }

    doubleConfParameter("V_rf",V_rf);
    quadFactor = - (N_rods-1.0)/(2*std::pow(rho,3.0)) * std::pow( N_rods * V_rf / omega_rf  ,2.0);

    doubleConfParameter("space_charge_factor", spaceChargeFactor);
    doubleConfParameter("collision_gas_mass_amu", collisionGasMassAmu);
    doubleConfParameter("start_z_length_mm", startZLength);
    startZLength = startZLength / 1000.0;

    doubleConfParameter("restart_z_length_mm", reStartZLength);
    reStartZLength = reStartZLength / 1000.0;

    doubleConfParameter("max_z_length_mm", maxZLength);
    maxZLength = maxZLength / 1000.0;

    doubleConfParameter("entrance_aperture_mm", entranceAperture);
    entranceAperture = entranceAperture / 1000.0;
    rnd_xy = std::uniform_real_distribution<double>(-entranceAperture,entranceAperture);

    doubleConfParameter("dist_entrance_mm", distEntrance);
    distEntrance = distEntrance / 1000.0;

    doubleConfParameter("dist_quad_mm", distQuadLength);
    distQuadLength = distQuadLength / 1000.0;

    doubleConfParameter("z_gradient_entrance_V", z_gradient_entrance);
    if (distEntrance>0){
        z_gradient_entrance = z_gradient_entrance / distEntrance;
    }

    doubleConfParameter("z_gradient_quad_V", z_gradient_quad);
    if (distQuadLength>0){
        z_gradient_quad = z_gradient_quad / distQuadLength;
    }

    doubleConfParameter("background_gas_z_velocity", backgroundGasZVelocity);
    backgroundGasVelocity = Core::Vector(0.0,0.0,backgroundGasZVelocity);
    doubleConfParameter("background_gas_base_pressure", backgroundBasePressure);
    doubleConfParameter("background_gas_entrance_pressure", backgroundEntrancePressure);
    doubleConfParameter("linear_pressure_drop_length_mm", linearPressureDropLength);
    linearPressureDropLength = linearPressureDropLength / 1000.0;



    std::vector<Core::Particle*>particles= std::vector<Core::Particle*>();

    for (int i=0; i<nIons.size(); i++){
        int nParticles = nIons[i];
        double mass = ionMasses[i];
        Core::Particle* ions= ParticleSimulation::util::getRandomIonsInBox(
                nParticles,1.0,
                Core::Vector(-entranceAperture/2.0,-entranceAperture/2.0,0),
                Core::Vector(entranceAperture,entranceAperture,startZLength)
        );
        for (int j=0; j<nParticles; j++){
            ions[j].setMassAMU(mass);
            particles.push_back(&ions[j]);
        }
    }

    CollisionModel::HardSphereModel hsModel = CollisionModel::HardSphereModel(0.0,298,collisionGasMassAmu);
    hsModel.setPressureFunction(backgroundGasPressureFunction);
    hsModel.setVelocityFunction(backgroundGasVelocityFunction);

    std::function<void(Core::Vector& newPartPos,Core::Particle* particle, int particleIndex, SpaceCharge::FieldCalculator& scFieldCalculator, double time,int timestep)> otherActionsFct;

    if (simType==2){
        otherActionsFct = otherActionsFunctionQuad_2dSim;
    }else if(simType==3){
        otherActionsFct = otherActionsFunctionQuad_3dSim;
    }

    AppUtils::Stopwatch stopWatch;
    stopWatch.start();
    ParticleSimulation::verletIntegrationSimulation(
            particles,dt,timeSteps,accelerationFunctionQuad, timestepWriteFunction,otherActionsFct,hsModel,projectName);

    stopWatch.stop();
    std::cout << particles[0]->getLocation()<<std::endl;
    std::cout << "elapsed wall time:"<< stopWatch.elapsedSecondsWall()<<std::endl;
    std::cout << "elapsed cpu time:"<< stopWatch.elapsedSecondsCPU()<<std::endl;



    return 0;
}



double spaceChargeFactor = 0.0;
double collisionGasMassAmu = 0.0;
double distEntrance = 0.0;//10/1000.0; //distance between entrance and quad
double entranceAperture = 0.0;
double distQuadLength = 0.0; //50/1000.0;
double z_gradient_entrance = 0.0; //5.0 / distEntranceQuad; // potential gradient between entrance and quad
double z_gradient_quad = 0.0; //1.0 / distQuadLength;  // potential gradient over the quad

double backgroundGasZVelocity = 0.0;
double backgroundBasePressure = 0.0; //0.7;
double backgroundEntrancePressure = 10; //30 from DSMC
double linearPressureDropLength = 60;

double startZLength = 0.0;
double reStartZLength = 0.0;
double maxZLength = 0.0;

void timestepWriteFunction(SpaceCharge::FieldCalculator& scFieldCalculator, std::vector<int> ionKeys, double time, int timestep){
    if (timestep % trajectoryWriteInterval ==0){
        
        std::cout<<"ts:"<<timestep<<" time:"<<time<<std::endl;
        char buffer [100];
        int cx;
        cx = snprintf ( buffer, 100, "%06d", timestep );
        //Core::SimpleVTKwriter vtkWriter = Core::SimpleVTKwriter(projectName+"_"+std::string(buffer));
        //vtkWriter.write(tree,false);
        jsonWriter->writeTimestep(tree,ionKeys,time,false);
        //jsonWriter->writeTimestep(tree,ionKeys, time,false);
    }
    if (timestep == timeSteps -1){
        
        jsonWriter->writeTimestep(tree,ionKeys, time,true);
        //jsonWriter->writeTimestep(tree,ionKeys, time,true);
        
        jsonWriter->writeSplatTimes(tree,ionKeys);
        jsonWriter->writeIonMasses(tree, ionKeys);
        
        std::cout<<"finished ts:"<<timestep<<" time:"<<time<<std::endl;
    }
}

Core::Vector accelerationFunctionQuad(Core::Particle* particle, int particleIndex, SpaceCharge::FieldCalculator& scFieldCalculator, double time,int timestep){
    //z is the long quad axis
    Core::Vector pos = particle->getLocation();
    double particleCharge =particle->getCharge();
    double r_pos = std::sqrt( pos.x()*pos.x() + pos.y()*pos.y() );
    double E_quad = quadFactor*( particleCharge / particle->getMass()) * std::pow(r_pos / rho, quadExp);
    
    double z_gradient = 0.0;
    if (pos.z() < distEntrance){
        z_gradient = z_gradient_entrance;
    }else{
        z_gradient = z_gradient_quad;
    }
   
    //

    Core::Vector spaceChargeForce = scFieldCalculator.computeEFieldFromSpaceCharge(*particle)* spaceChargeFactor;

    
    Core::Vector result = (Core::Vector(
                                          E_quad * pos.x(),
                                          E_quad * pos.y(),
                                          z_gradient)  + spaceChargeForce) * particleCharge / particle->getMass();
    //std::cout << result << "  "<<pos<<" "<<particleCharge / particle->getMass()<< std::endl;
    
    return(result);
}

Core::Vector backgroundGasVelocity= Core::Vector(0.0,0.0,0.0);
Core::Vector backgroundGasVelocityFunction(Core::Vector& location){
    return backgroundGasVelocity;
}

double backgroundGasPressureFunction(Core::Vector& location){
    double zPos = location.z();
    if (zPos > linearPressureDropLength){
        return backgroundBasePressure;
    }
    else if(zPos < 0.0){
        return backgroundEntrancePressure;
    }
    else {
        return ( (1 - (zPos / linearPressureDropLength)) * (backgroundEntrancePressure-backgroundBasePressure) + backgroundBasePressure );
    }
}

std::uniform_real_distribution<double> rnd_xy(0.0,0.0);
std::uniform_real_distribution<double> rnd_z(0,1.0);

void otherActionsFunctionQuad_3dSim(Core::Vector& newPartPos,Core::Particle* particle, int particleIndex, SpaceCharge::FieldCalculator& scFieldCalculator, double time,int timestep){
    double r_pos = std::sqrt( newPartPos.x()*newPartPos.x() + newPartPos.y()*newPartPos.y() );
    
/*    if (newPartPos.z() > maxZLength){
        newPartPos.z(0.0);
    }
 */
    if (r_pos > rho || newPartPos.z() > maxZLength ){
        double xNew = rnd_xy(BTree::globalRNG);
        double yNew = rnd_xy(BTree::globalRNG);
        
        newPartPos.set(
                xNew,//rnd_xy(Core::globalRNG),
                yNew,//rnd_xy(Core::globalRNG),
                rnd_z (BTree::globalRNG)*reStartZLength);
    }
}

void otherActionsFunctionQuad_2dSim(Core::Vector& newPartPos,Core::Particle* particle, int particleIndex, SpaceCharge::FieldCalculator& scFieldCalculator, double time,int timestep){
    newPartPos.z(startZLength);
}

