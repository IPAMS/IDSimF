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

#include "RS_Simulation.hpp"
#include "RS_SimulationConfiguration.hpp"
#include "RS_ConfigFileParser.hpp"
#include "RS_ConcentrationFileWriter.hpp"
#include "PSim_trajectoryExplorerJSONwriter.hpp"
#include "PSim_sampledWaveform.hpp"
#include "PSim_scalar_writer.hpp"
#include "PSim_util.hpp"
#include "appUtils_simulationConfiguration.hpp"
#include "appUtils_logging.hpp"
#include "appUtils_stopwatch.hpp"
#include "appUtils_signalHandler.hpp"
#include "json.h"
#include <iostream>
#include <cmath>

const std::string key_ChemicalIndex = "keyChemicalIndex";
const double standardPressure_Pa = 102300; //101325
enum CVMode {STATIC_CV, AUTO_CV, MODULATED_CV, MODULATED_AUTO_CV};
enum SVMode {BI_SIN, SQUARE, CLIPPED_SIN};

int main(int argc, const char * argv[]) {

    try{
        // open configuration, parse configuration file =========================================
        if (argc<=2) {
            std::cout << "Run abort: No run configuration or project name given." << std::endl;
            return EXIT_FAILURE;
        }
        std::string projectName = argv[2];
        std::cout << projectName << std::endl;
        auto logger = AppUtils::createLogger(projectName + ".log");

        std::string confFileName = argv[1];
        AppUtils::SimulationConfiguration simConf(confFileName, logger);

        std::vector<unsigned int> nParticles = simConf.unsignedIntVectorParameter("n_particles");
        int nSteps = simConf.intParameter("sim_time_steps");
        int nStepsPerOscillation = simConf.intParameter("sim_time_steps_per_sv_oscillation");
        int concentrationWriteInterval = simConf.intParameter("concentrations_write_interval");
        int trajectoryWriteInterval = simConf.intParameter("trajectory_write_interval");

        //geometric parameters:
        double startWidthX_m = simConf.doubleParameter("start_width_x_mm")/1000.0;
        double startWidthY_m = simConf.doubleParameter("start_width_y_mm")/1000.0;
        double startWidthZ_m = simConf.doubleParameter("start_width_z_mm")/1000.0;


        //Define background temperature
        double backgroundTemperature_K = simConf.doubleParameter("background_temperature_K");
        double backgroundPressure_Pa = simConf.doubleParameter("background_pressure_Pa");

        //field parameters:
        std::string cvModeStr = simConf.stringParameter("cv_mode");
        CVMode cvMode;
        double meanZPos = 0.0; //variable used for automatic CV correction
        double cvRelaxationParameter = 0.0;

        if (cvModeStr == "static"){
            cvMode = STATIC_CV;
        }
        else if (cvModeStr == "auto"){
            cvMode = AUTO_CV;
            cvRelaxationParameter = simConf.doubleParameter("cv_relaxation_parameter");
        }
        else if (cvModeStr == "modulated"){
            cvMode = MODULATED_CV;
        }
        else if (cvModeStr == "modulated_auto"){
            cvMode = MODULATED_AUTO_CV;
        }
        else{
            throw std::invalid_argument("wrong configuration value: cv_mode");
        }


        SVMode svMode;
        std::string svModeStr = simConf.stringParameter("sv_mode");
        if (svModeStr == "bi_sin"){
            svMode = BI_SIN;
        }
        else if (svModeStr == "square"){
            svMode = SQUARE;
        }
        else if (svModeStr == "clipped_sin"){
            svMode = CLIPPED_SIN;
        }
        else{
            throw std::invalid_argument("wrong configuration value: sv_mode");
        }

        double fieldSVSetpoint_VPerM = simConf.doubleParameter("sv_Vmm-1") * 1000.0;
        double fieldCVSetpoint_VPerM = simConf.doubleParameter("cv_Vmm-1") * 1000.0;
        double fieldFrequency = simConf.doubleParameter("sv_frequency_s-1");
        double fieldWavePeriod = 1.0/fieldFrequency;
        double fieldMagnitude = 0; //actual field magnitude
        double dt_s = fieldWavePeriod / nStepsPerOscillation;

        // ======================================================================================

        //read and prepare chemical configuration ===============================================
        RS::ConfigFileParser parser = RS::ConfigFileParser();
        std::string rsConfFileName = simConf.pathRelativeToConfFile(simConf.stringParameter("reaction_configuration"));
        RS::Simulation rsSim = RS::Simulation(parser.parseFile(rsConfFileName));
        RS::SimulationConfiguration* rsSimConf = rsSim.simulationConfiguration();
        //prepare a map for retrieval of the substance index:
        std::map<RS::Substance*,int> substanceIndices;
        std::vector<RS::Substance*> discreteSubstances = rsSimConf->getAllDiscreteSubstances();
        std::vector<double> ionMobility; // = simConf.doubleVectorParameter("ion_mobility");
        for (std::size_t i=0; i<discreteSubstances.size(); i++){
            substanceIndices.insert(std::pair<RS::Substance*,int>(discreteSubstances[i], i));
            ionMobility.push_back(discreteSubstances[i]->mobility());
        }


        // prepare file writer  =================================================================
        RS::ConcentrationFileWriter resultFilewriter(projectName+"_conc.csv");

        auto jsonWriter = std::make_unique<ParticleSimulation::TrajectoryExplorerJSONwriter>(projectName+ "_trajectories.json");
        jsonWriter->setScales(1000,1);

        // read particle configuration ==========================================================
        unsigned int ionsInactive = 0;
        unsigned int nAllParticles = 0;
        for (const auto ni: nParticles){
            nAllParticles += ni;
        }

        std::unique_ptr<ParticleSimulation::Scalar_writer> cvFieldWriter;
        int cvHighResLogPeriod = 0;
        if (simConf.isParameter("log_cv_field_period")){
            cvHighResLogPeriod = simConf.intParameter("log_cv_field_period");
        }
        if (cvMode == AUTO_CV || cvHighResLogPeriod > 0){
            cvFieldWriter = std::make_unique<ParticleSimulation::Scalar_writer>(projectName+ "_cv.csv");
        }

        std::unique_ptr<ParticleSimulation::Scalar_writer> voltageWriter;
        voltageWriter = std::make_unique<ParticleSimulation::Scalar_writer>(projectName+ "_voltages.csv");

        // init simulation  =====================================================================

        // create and add simulation particles:
        unsigned int nParticlesTotal = 0;
        std::vector<uniqueReactivePartPtr>particles;
        std::vector<BTree::Particle*>particlesPtrs;
        std::vector<std::vector<double>> trajectoryAdditionalParams;

        Core::Vector initCorner(0,-startWidthY_m/2.0,-startWidthZ_m/2.0);
        Core::Vector initBoxSize(startWidthX_m,startWidthY_m,startWidthZ_m);

        for (std::size_t i=0; i<nParticles.size(); ++i) {
            RS::Substance *subst = rsSimConf->substance(i);
            std::vector<Core::Vector> initialPositions =
                    ParticleSimulation::util::getRandomPositionsInBox(nParticles[i],initCorner,initBoxSize);
            for (unsigned int k = 0; k < nParticles[i]; k++) {
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

        resultFilewriter.initFile(rsSimConf);
        // ======================================================================================


        // define trajectory integration parameters / functions =================================
        std::function<double(double fieldAmplitude_VPerM, double time)> SVFieldFct;

        if (svMode == BI_SIN){
            double field_h = 2.0;
            double field_F = 2.0;
            double field_W = (1/fieldWavePeriod) * 2 * M_PI;
            auto  fieldFctBisinusoidal=
                    [field_F, field_W, field_h]
                    (double svAmplitude_VPerM, double time) -> double{

                        //double particleCharge = particle->getCharge();
                        double voltageSVgp = svAmplitude_VPerM * 0.6667; // V/m (1V/m peak to peak is 0.6667V/m ground to peak)
                        double voltageSVt = (field_F * sin(field_W * time)
                                + sin(field_h * field_W * time - 0.5 * M_PI))
                                * voltageSVgp / (field_F + 1);

                        return voltageSVt;
                    };
            SVFieldFct = fieldFctBisinusoidal;
        }
        else if (svMode == SQUARE){
            double thirdOfWavePeriod = fieldWavePeriod / 3.0;
            auto  fieldFctSquare=
                    [fieldWavePeriod, thirdOfWavePeriod]
                    (double svAmplitude_VPerM, double time) -> double{

                        double timeInPeriod = std::fmod(time, fieldWavePeriod);

                        double voltageSVgp_highField = svAmplitude_VPerM * 0.666667; // V/m (1V/m peak to peak is 0.6667V/m ground to peak)
                        double voltageSVgp_lowField = -voltageSVgp_highField * 0.5; // low field is 1/2 of high field
                        if (timeInPeriod < thirdOfWavePeriod){
                            return voltageSVgp_highField;
                        }
                        else {
                            return voltageSVgp_lowField;
                        }
                    };
            SVFieldFct = fieldFctSquare;
        }
        else if (svMode == CLIPPED_SIN){
            double h = 3.0/2.0* M_PI  - 1.0;
            double t_sin = M_PI / (2*h + 1);
            double f_low = -(2.0* t_sin) / (M_PI-2.0*t_sin);

            auto  fieldFctClippedSin=
                    [f_low, t_sin, fieldWavePeriod]
                    (double svAmplitude_VPerM, double time) -> double{

                        double normalizedTimeInPeriod = std::fmod(time, fieldWavePeriod) / fieldWavePeriod;
                        if (normalizedTimeInPeriod < t_sin){
                            return  ((M_PI * sin(M_PI* normalizedTimeInPeriod / t_sin) - 2*t_sin) / (M_PI-2*t_sin) * svAmplitude_VPerM);
                        }
                        else {
                            return  (f_low*svAmplitude_VPerM);
                        }
                    };
            SVFieldFct = fieldFctClippedSin;
        }

        std::function<double(double cvAmplitude_VPerM, double time)> CVFieldFct;
        if (cvMode == STATIC_CV || cvMode == AUTO_CV){
            //non modulated CV:
            auto nonModulatedCV =
                [] (double cvAmplitude_VPerM, double) -> double{
                    return cvAmplitude_VPerM;
                };
            CVFieldFct = nonModulatedCV;
        }
        else {
            // modulated CV, read CV waveform and phase shift:
            std::string cvWaveformFileName = simConf.pathRelativeToConfFile(simConf.stringParameter("cv_waveform"));
            auto cvWaveForm = ParticleSimulation::SampledWaveform(cvWaveformFileName);

            double cvPhaseShift = simConf.doubleParameter("cv_phase_shift");
            auto modulatedCV =
                [cvWaveForm, cvPhaseShift, fieldWavePeriod] (double cvAmplitude_VPerM, double time) -> double{
                    double period = std::fmod(time, fieldWavePeriod) / fieldWavePeriod;
                    double shiftedPeriod = std::fmod(period + cvPhaseShift, 1.0);
                    return cvWaveForm.getInterpolatedValue(shiftedPeriod) * cvAmplitude_VPerM;
                };
            CVFieldFct = modulatedCV;
        }

        ParticleSimulation::partAttribTransformFctType additionalParameterTransformFct =
                [=](BTree::Particle* particle) -> std::vector<double> {
                    std::vector<double> result = {particle->getFloatAttribute(key_ChemicalIndex)};
                    return result;
                };

        auto timestepWriteFct =
                [&jsonWriter, &voltageWriter, &additionalParameterTransformFct, trajectoryWriteInterval,
                        &rsSim, &resultFilewriter, concentrationWriteInterval, &fieldMagnitude, &logger]
                        (std::vector<BTree::Particle *> &particles, double time, int timestep, bool lastTimestep){

            if (timestep % concentrationWriteInterval ==0) {
                resultFilewriter.writeTimestep(rsSim);
                voltageWriter->writeTimestep(fieldMagnitude, time);
            }
            if (lastTimestep) {
                jsonWriter->writeTimestep(
                        particles, additionalParameterTransformFct, time, true);

                jsonWriter->writeSplatTimes(particles);
                jsonWriter->writeIonMasses(particles);
                logger->info("finished ts:{} time:{:.2e}", timestep, time);
            }

            else if (timestep % trajectoryWriteInterval ==0){
                rsSim.logConcentrations(logger);
                jsonWriter->writeTimestep(
                        particles, additionalParameterTransformFct, time, false);
            }
        };


        // simulate   ===========================================================================
        AppUtils::SignalHandler::registerSignalHandler();
        AppUtils::Stopwatch stopWatch;
        stopWatch.start();

        reactionConditions.temperature = backgroundTemperature_K;
        double reducedPressure = standardPressure_Pa / backgroundPressure_Pa;
        for (int step=0; step<nSteps; step++) {

            double cvFieldNow_VPerM = CVFieldFct(fieldCVSetpoint_VPerM, rsSim.simulationTime());
            double svFieldNow_VPerM = SVFieldFct(fieldSVSetpoint_VPerM, rsSim.simulationTime());
            fieldMagnitude = svFieldNow_VPerM + cvFieldNow_VPerM;
            reactionConditions.electricField = fieldMagnitude;

            for (unsigned int i = 0; i < nParticlesTotal; i++) {
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
            if (cvHighResLogPeriod > 0 && step % cvHighResLogPeriod == 0){
                cvFieldWriter->writeTimestep(cvFieldNow_VPerM, rsSim.simulationTime());
            }
            rsSim.advanceTimestep(dt_s);
            timestepWriteFct(particlesPtrs, rsSim.simulationTime(), step, false);

            //autocorrect compensation voltage, to minimize z drift (once for every single SV oscillation):
            if (cvMode == AUTO_CV && step % nStepsPerOscillation == 0) {
                //calculate current mean z-position:
                double buf = 0.0;
                for (unsigned int i = 0; i < nParticlesTotal; i++) {
                    buf += particles[i]->getLocation().z();
                }
                double currentMeanZPos = buf / nParticlesTotal;

                //update cv value:
                double diffMeanZPos = meanZPos - currentMeanZPos;
                fieldCVSetpoint_VPerM = fieldCVSetpoint_VPerM + diffMeanZPos * cvRelaxationParameter;
                cvFieldWriter->writeTimestep(std::vector<double>{fieldCVSetpoint_VPerM, currentMeanZPos}, rsSim.simulationTime());
                meanZPos = currentMeanZPos;
                logger->info("CV corrected ts:{} time:{:.2e} new CV:{}", step, rsSim.simulationTime(), fieldCVSetpoint_VPerM);
            }

            if (ionsInactive>=nAllParticles || AppUtils::SignalHandler::isTerminationSignaled()){
                timestepWriteFct(particlesPtrs, rsSim.simulationTime(), step, true);
                break;
            }
        }
        resultFilewriter.closeFile();

        stopWatch.stop();

        logger->info("total reaction events: {} ill events: {}", rsSim.totalReactionEvents(), rsSim.illEvents());
        logger->info("ill fraction: {}", rsSim.illEvents() / (double) rsSim.totalReactionEvents());
        logger->info("CPU time: {} s", stopWatch.elapsedSecondsCPU());
        logger->info("Finished in {} seconds (wall clock time)",stopWatch.elapsedSecondsWall());
        // ======================================================================================

        return 0;
    }
    catch(const std::invalid_argument& ia){
        std::cout << ia.what() << std::endl;
        return EXIT_FAILURE;
    }
}
