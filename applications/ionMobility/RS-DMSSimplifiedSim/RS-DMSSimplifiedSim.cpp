/**************************
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
 BT-RS-DMSSimplifiedSim.cpp

 Idealized plane electrode type differential ion mobility spectrometry (DMS):
 Simplified / idealized simulation of compensation voltage with chemically reactive
 particle ensemble

 ****************************/

#include <iostream>
#include <cmath>
#include <ctime>
#include "json.h"
#include "appUtils_parameterParsing.hpp"
#include "RS_Simulation.hpp"
#include "RS_SimulationConfiguration.hpp"
#include "RS_ConfigFileParser.hpp"
#include "RS_ConcentrationFileWriter.hpp"
#include "PSim_trajectoryExplorerJSONwriter.hpp"
#include "PSim_scalar_writer.hpp"
#include "PSim_util.hpp"

const std::string key_ChemicalIndex = "keyChemicalIndex";
const double standardPressure_Pa = 102300; //101325
enum CVMode {STATIC_CV, AUTO_CV};

int main(int argc, const char * argv[]) {

    // open configuration, parse configuration file =========================================
    if (argc <2){
        std::cout << "no configuration file given"<<std::endl;
        return(1);
    }

    std::string confFileName = argv[1];
    Json::Value confRoot = readConfigurationJson(confFileName);

    std::vector<int> nParticles = intVectorConfParameter("n_particles",confRoot);
    int nSteps = intConfParameter("sim_time_steps",confRoot);
    int nStepsPerOscillation = intConfParameter("sim_time_steps_per_sv_oscillation",confRoot);
    int concentrationWriteInterval = intConfParameter("concentrations_write_interval",confRoot);
    int trajectoryWriteInterval = intConfParameter("trajectory_write_interval",confRoot);

    //geometric parameters:
    double startWidthX_m = doubleConfParameter("start_width_x_mm",confRoot)/1000.0;
    double startWidthY_m = doubleConfParameter("start_width_y_mm",confRoot)/1000.0;
    double startWidthZ_m = doubleConfParameter("start_width_z_mm",confRoot)/1000.0;


    //Define background temperature
    double backgroundTemperature_K = doubleConfParameter("background_temperature_K",confRoot);
    double backgroundPressure_Pa = doubleConfParameter("background_pressure_Pa",confRoot);

    //field parameters:
    std::string cvModeStr = stringConfParameter("cv_mode",confRoot);
    CVMode cvMode;
    double meanZPos = 0.0; //variable used for automatic CV correction
    double cvRelaxationParameter;
    if (cvModeStr == "static"){
        cvMode = STATIC_CV;
    }
    else if (cvModeStr == "auto"){
        cvMode = AUTO_CV;
        cvRelaxationParameter = doubleConfParameter("cv_relaxation_parameter",confRoot);
    }
    else{
        throw std::invalid_argument("wrong configuration value: cv_mode");
    }

    double fieldSV_VPerM = doubleConfParameter("sv_Vmm-1",confRoot) * 1000.0;
    double fieldCV_VPerM = doubleConfParameter("cv_Vmm-1",confRoot) * 1000.0;
    double fieldFrequency = doubleConfParameter("sv_frequency_s-1",confRoot);
    double fieldWavePeriod = 1.0/fieldFrequency;
    double field_h = 2.0;
    double field_F = 2.0;
    double field_W = (1/fieldWavePeriod) * 2 * M_PI;
    double fieldMagnitude = 0; //actual field magnitude

    double dt_s = fieldWavePeriod / nStepsPerOscillation;

    std::string projectFilename;
    if (argc == 3){
        projectFilename = argv[2];
    }
    else {
        std::stringstream ss;
        ss << argv[1] << "_result.txt";
        projectFilename = ss.str();
    }
    // ======================================================================================


    //read and prepare chemical configuration ===============================================
    RS::ConfigFileParser parser = RS::ConfigFileParser();
    std::string rsConfFileName = pathRelativeToConfFile(
                                    confFileName,
                                    stringConfParameter("reaction_configuration",confRoot));
    RS::Simulation rsSim = RS::Simulation(parser.parseFile(rsConfFileName));
    RS::SimulationConfiguration* simConf = rsSim.simulationConfiguration();
    //prepare a map for retrieval of the substance index:
    std::map<RS::Substance*,int> substanceIndices;
    std::vector<RS::Substance*> discreteSubstances = simConf->getAllDiscreteSubstances();
    std::vector<double> ionMobility; // = doubleVectorConfParameter("ion_mobility",confRoot);
    for (int i=0; i<discreteSubstances.size(); i++){
        substanceIndices.insert(std::pair<RS::Substance*,int>(discreteSubstances[i], i));
        ionMobility.push_back(discreteSubstances[i]->mobility());
    }

    // prepare file writer  =================================================================
    RS::ConcentrationFileWriter resultFilewriter(projectFilename+"_conc.csv");

    auto jsonWriter = std::make_unique<ParticleSimulation::TrajectoryExplorerJSONwriter>(projectFilename+ "_trajectories.json");
    jsonWriter->setScales(1000,1);

    // read particle configuration ==========================================================
    int ionsInactive = 0;
    int nAllParticles = 0;
    for (const auto ni: nParticles){
        nAllParticles += ni;
    }

    std::unique_ptr<ParticleSimulation::Scalar_writer> cvFieldWriter;
    if (cvMode == AUTO_CV){
        cvFieldWriter = std::make_unique<ParticleSimulation::Scalar_writer>(projectFilename+ "_cv.csv");
    }

    std::unique_ptr<ParticleSimulation::Scalar_writer> voltageWriter;
    voltageWriter = std::make_unique<ParticleSimulation::Scalar_writer>(projectFilename+ "_voltages.csv");


    // init simulation  =====================================================================

    // create and add simulation particles:
    int nParticlesTotal = 0;
    std::vector<uniqueReactivePartPtr>particles;
    std::vector<BTree::Particle*>particlesPtrs;
    std::vector<std::vector<double>> trajectoryAdditionalParams;

    Core::Vector initCorner(0,-startWidthY_m/2.0,-startWidthZ_m/2.0);
    Core::Vector initBoxSize(startWidthX_m,startWidthY_m,startWidthZ_m);

    for (int i=0; i<nParticles.size();i++) {
        RS::Substance *subst = simConf->substance(i);
        std::vector<Core::Vector> initialPositions =
                ParticleSimulation::util::getRandomPositionsInBox(nParticles[i],initCorner,initBoxSize);
        for (int k = 0; k < nParticles[i]; k++) {
            uniqueReactivePartPtr particle = std::make_unique<RS::ReactiveParticle>(subst);

            // init position and initial chemical species of the particle:
            particle->setLocation(initialPositions[k]);
            int substIndex = substanceIndices.at(particle->getSpecies());
            particle->setFloatAttribute(key_ChemicalIndex, substIndex);

            particlesPtrs.push_back(particle.get());
            rsSim.addParticle(particle.get(), nParticlesTotal);
            particles.push_back(std::move(particle));
            trajectoryAdditionalParams.emplace_back(std::vector<double>(1));
            nParticlesTotal++;
        }
    }

    RS::ReactionConditions reactionConditions = RS::ReactionConditions();
    reactionConditions.temperature = 0.0;//backgroundTemperature_K;
    reactionConditions.pressure = backgroundPressure_Pa;
    reactionConditions.electricField = 0.0;
    reactionConditions.totalReactionEnergy = 0.0;

    resultFilewriter.initFile(simConf);
    // ======================================================================================


    // define trajectory integration parameters / functions =================================

    auto fieldFct =
            [&fieldSV_VPerM, &fieldCV_VPerM, field_F, field_W, field_h]
            (double time) -> double{

        //double particleCharge = particle->getCharge();
        double voltageSVgp = fieldSV_VPerM * 0.6667; // V/m (1V/m peak to peak is 0.6667V/m ground to peak)
        double voltageSVt = fieldCV_VPerM + (field_F * sin(field_W * time)
                                    + sin(field_h * field_W * time - 0.5 * M_PI))
                                    * voltageSVgp / (field_F + 1);

        //fieldMagnitude = voltageSVt; //// / electrodeDistance_m;

        //Core::Vector fieldForce(0, 0, fieldMagnitude * particleCharge);

        return voltageSVt;
    };


    ParticleSimulation::additionalPartParamFctType additionalParameterTransformFct =
            [=](BTree::Particle* particle) -> std::vector<double> {
                std::vector<double> result = {particle->getFloatAttribute(key_ChemicalIndex)};
                return result;
            };


    auto timestepWriteFct =
            [&jsonWriter, &voltageWriter, &additionalParameterTransformFct, trajectoryWriteInterval,
                    &rsSim, &resultFilewriter, concentrationWriteInterval, &fieldMagnitude]
                    (std::vector<BTree::Particle *> &particles, double time, int timestep, bool lastTimestep){

        if (timestep % concentrationWriteInterval ==0) {
            rsSim.printConcentrations();
            resultFilewriter.writeTimestep(rsSim);
            voltageWriter->writeTimestep(fieldMagnitude,time);
        }
        if (lastTimestep) {
            jsonWriter->writeTimestep(
                    particles, additionalParameterTransformFct, time, true);

            jsonWriter->writeSplatTimes(particles);
            jsonWriter->writeIonMasses(particles);

            std::cout << "finished ts:" << timestep << " time:" << time << std::endl;
        }

        else if (timestep % trajectoryWriteInterval ==0){
            jsonWriter->writeTimestep(
                    particles, additionalParameterTransformFct, time, false);
        }
    };


    // simulate   ===========================================================================
    clock_t begin = std::clock();

    reactionConditions.temperature = backgroundTemperature_K;
    double reducedPressure = standardPressure_Pa / backgroundPressure_Pa;
    for (int step=0; step<nSteps; step++) {

        fieldMagnitude = fieldFct(rsSim.simulationTime());
        reactionConditions.electricField = fieldMagnitude;

        for (int i = 0; i < nParticlesTotal; i++) {
            bool reacted = rsSim.react(i, reactionConditions, dt_s);

            if (reacted){
                //we had an reaction event: Update the chemical species for the trajectory
                int substIndex = substanceIndices.at(particles[i]->getSpecies());
                particles[i]->setFloatAttribute(key_ChemicalIndex, substIndex);
            }
            // move particle in z axis according to particle mobility:
            double particleShift = particles[i]->getMobility() * reducedPressure * fieldMagnitude * dt_s;
            Core::Vector particleLocation = particles[i]->getLocation();
            particleLocation.z( particleLocation.z() + particleShift);
            particles[i]->setLocation(particleLocation);
        }
        rsSim.advanceTimestep(dt_s);
        timestepWriteFct(particlesPtrs, rsSim.simulationTime(), step, false);

        //autocorrect compensation voltage, to minimize z drift (once for every single SV oscillation):
        if (cvMode == AUTO_CV && step % nStepsPerOscillation == 0) {
            //calculate current mean z-position:
            double buf = 0.0;
            for (int i = 0; i < nParticlesTotal; i++) {
                buf += particles[i]->getLocation().z();
            }
            double currentMeanZPos = buf / nParticlesTotal;

            //update cv value:
            double diffMeanZPos = meanZPos - currentMeanZPos;
            fieldCV_VPerM = fieldCV_VPerM + diffMeanZPos * cvRelaxationParameter;
            cvFieldWriter->writeTimestep(std::vector<double>{fieldCV_VPerM,currentMeanZPos},rsSim.simulationTime());
            meanZPos = currentMeanZPos;
            std::cout << "CV corrected ts:"<< step << " time:"<<rsSim.simulationTime()<< " new val:"<<fieldCV_VPerM<<std::endl;
        }

        if (ionsInactive >= nAllParticles){
            break;
        }
    }
    timestepWriteFct(particlesPtrs, rsSim.simulationTime(), nSteps, true);

    clock_t end = std::clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "reaction events:" << rsSim.totalReactionEvents() << " ill events:" << rsSim.illEvents() << std::endl;
    std::cout << "ill fraction: " << rsSim.illEvents() / (double) rsSim.totalReactionEvents() << std::endl;
    std::cout << "elapsed seconds "<< elapsed_secs<<std::endl;

    // ======================================================================================


    return 0;
}
