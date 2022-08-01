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
#include "FileIO_trajectoryHDF5Writer.hpp"
#include "FileIO_scalar_writer.hpp"
#include "PSim_util.hpp"
#include "PSim_sampledFunction.hpp"
#include "appUtils_simulationConfiguration.hpp"
#include "appUtils_logging.hpp"
#include "appUtils_stopwatch.hpp"
#include "appUtils_signalHandler.hpp"
#include "appUtils_commandlineParser.hpp"
#include "dmsSim_dmsFields.hpp"
#include "json.h"
#include <iostream>
#include <cmath>

const std::string key_ChemicalIndex = "keyChemicalIndex";
const double standardPressure_Pa = 102300; //101325
enum MobilityScalingMode {STATIC_MOBILITY, FIELD_SCALING_FUNCTION};

int main(int argc, const char * argv[]) {

    try{
        AppUtils::CommandlineParser cmdLineParser(argc, argv, "DMS simplified", "Simplified DMS Simulation", true);
        std::string projectName = cmdLineParser.resultName();
        AppUtils::logger_ptr logger = cmdLineParser.logger();

        std::string confFileName = cmdLineParser.confFileName();
        AppUtils::simConf_ptr simConf = cmdLineParser.simulationConfiguration();

        std::vector<unsigned int> nParticles = simConf->unsignedIntVectorParameter("n_particles");
        int nSteps = simConf->intParameter("sim_time_steps");
        int nStepsPerOscillation = simConf->intParameter("sim_time_steps_per_sv_oscillation");
        int concentrationWriteInterval = simConf->intParameter("concentrations_write_interval");
        int trajectoryWriteInterval = simConf->intParameter("trajectory_write_interval");

        //geometric parameters:
        double startWidthX_m = simConf->doubleParameter("start_width_x_mm")/1000.0;
        double startWidthY_m = simConf->doubleParameter("start_width_y_mm")/1000.0;
        double startWidthZ_m = simConf->doubleParameter("start_width_z_mm")/1000.0;


        //Define background temperature
        double backgroundTemperature_K = simConf->doubleParameter("background_temperature_K");
        double backgroundPressure_Pa = simConf->doubleParameter("background_pressure_Pa");
        double reducedPressure = standardPressure_Pa / backgroundPressure_Pa;

        //field parameters:
        CVMode cvMode = parseCVModeConfiguration(simConf);
        double meanZPos = 0.0; //variable used for automatic CV correction
        double cvRelaxationParameter = 0.0;
        if (cvMode == AUTO_CV || cvMode == MODULATED_AUTO_CV){
            cvRelaxationParameter = simConf->doubleParameter("cv_relaxation_parameter");
        }

        SVMode svMode = parseSVModeConfiguration(simConf);

        double fieldSVSetpoint_VPerM = simConf->doubleParameter("sv_Vmm-1") * 1000.0;
        double fieldCVSetpoint_VPerM = simConf->doubleParameter("cv_Vmm-1") * 1000.0;
        double fieldFrequency = simConf->doubleParameter("sv_frequency_s-1");
        double fieldWavePeriod = 1.0/fieldFrequency;
        double fieldMagnitude = 0; //actual field magnitude
        double dt_s = fieldWavePeriod / nStepsPerOscillation;

        // ======================================================================================

        //read and prepare chemical configuration ===============================================
        RS::ConfigFileParser parser = RS::ConfigFileParser();
        std::string rsConfFileName = simConf->pathRelativeToConfFile(simConf->stringParameter("reaction_configuration"));
        RS::Simulation rsSim = RS::Simulation(parser.parseFile(rsConfFileName));
        RS::SimulationConfiguration* rsSimConf = rsSim.simulationConfiguration();
        //prepare a map for retrieval of the substance index:
        std::map<RS::Substance*, int> substanceIndices;
        std::vector<RS::Substance*> discreteSubstances = rsSimConf->getAllDiscreteSubstances();
        std::vector<double> ionMobility; // = simConf->doubleVectorParameter("ion_mobility");
        for (std::size_t i=0; i<discreteSubstances.size(); i++){
            substanceIndices.insert(std::pair<RS::Substance*,int>(discreteSubstances[i], i));
            ionMobility.push_back(discreteSubstances[i]->lowFieldMobility());
        }

        //read ion mobility scaling configuration ==============================================
        MobilityScalingMode mobilitScalingMode;
        std::unique_ptr<ParticleSimulation::SampledFunction> mobilityScalingFct= nullptr;
        if (simConf->isParameter("mobility_scaling_function")){
            std::string  scalingFunctionFilename = simConf->pathRelativeToConfFile(
                    simConf->stringParameter("mobility_scaling_function"));
            mobilityScalingFct = std::make_unique<ParticleSimulation::SampledFunction>(scalingFunctionFilename);
            mobilitScalingMode = FIELD_SCALING_FUNCTION;
            if(!mobilityScalingFct->good()){
                throw(std::invalid_argument("Mobility function data not good"));
            }
        }
        else{
            mobilitScalingMode = STATIC_MOBILITY;
        }


        // prepare file writer  =================================================================
        RS::ConcentrationFileWriter resultFilewriter(projectName+"_conc.csv");

        std::vector<std::string> auxParamNames = {"chemical id"};
        auto additionalParamTFct = [](Core::Particle* particle) -> std::vector<int> {
            std::vector<int> result = {
                    particle->getIntegerAttribute(key_ChemicalIndex)
            };
            return result;
        };

        std::string hdf5Filename = projectName+"_trajectories.hd5";
        FileIO::TrajectoryHDF5Writer trajectoryWriter(hdf5Filename);
        trajectoryWriter.setParticleAttributes(auxParamNames, additionalParamTFct);

        std::unique_ptr<FileIO::Scalar_writer> cvFieldWriter;
        int cvHighResLogPeriod = 0;
        if (simConf->isParameter("log_cv_field_period")){
            cvHighResLogPeriod = simConf->intParameter("log_cv_field_period");
        }
        if (cvMode == AUTO_CV || cvHighResLogPeriod > 0){
            cvFieldWriter = std::make_unique<FileIO::Scalar_writer>(projectName+ "_cv.csv");
        }

        std::unique_ptr<FileIO::Scalar_writer> voltageWriter;
        voltageWriter = std::make_unique<FileIO::Scalar_writer>(projectName+ "_voltages.csv");

        std::unique_ptr<FileIO::Scalar_writer> mobilityWriter;
        mobilityWriter = std::make_unique<FileIO::Scalar_writer>(projectName+ "_mobility.csv");


        // init simulation  =====================================================================

        // create and add simulation particles:
        // read particle configuration ==========================================================
        unsigned int ionsInactive = 0;
        unsigned int nParticlesTotal = 0;
        std::vector<uniqueReactivePartPtr>particles;
        std::vector<Core::Particle*>particlesPtrs;
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
                particle->setIntegerAttribute(key_ChemicalIndex, substIndex);

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
        SVFieldFctType SVFieldFct = createSVFieldFunction(svMode, fieldWavePeriod);
        CVFieldFctType CVFieldFct = createCVFieldFunction(cvMode, fieldWavePeriod, simConf);

        auto timestepWriteFct =
                [&trajectoryWriter, &voltageWriter, trajectoryWriteInterval, reducedPressure,
                        &rsSim, &resultFilewriter, &mobilityWriter, concentrationWriteInterval, &fieldMagnitude, &logger]
                        (std::vector<Core::Particle *> &particles, double time, int timestep, bool lastTimestep){

            if (timestep % concentrationWriteInterval ==0) {
                resultFilewriter.writeTimestep(rsSim);
                voltageWriter->writeTimestep(fieldMagnitude, time);

                // calculate average mobility
                // calculate summed mobility:
                double mobilitySum = 0.0;
                std::size_t nParticles = particles.size();

                #pragma omp parallel for reduction (+:mobilitySum) default(none) shared(particles) firstprivate(nParticles, reducedPressure)
                for (std::size_t i = 0; i<nParticles; ++i) {
                    mobilitySum += particles[i]->getMobility() * reducedPressure;
                }

                double averageMobility = mobilitySum/nParticles;
                mobilityWriter->writeTimestep(averageMobility, time);

            }
            if (lastTimestep) {
                trajectoryWriter.writeTimestep(particles, time);
                trajectoryWriter.writeSplatTimes(particles);
                trajectoryWriter.finalizeTrajectory();
                trajectoryWriter.writeTimestep(particles, time);
            }

            else if (timestep % trajectoryWriteInterval ==0){
                rsSim.logConcentrations(logger);
                trajectoryWriter.writeTimestep(particles, time);
            }
        };

        auto particlesReactedFct = [substanceIndices](RS::ReactiveParticle* particle){
            //we had an reaction event: Update the chemical species for the trajectory
            int substIndex = substanceIndices.at(particle->getSpecies());
            particle->setIntegerAttribute(key_ChemicalIndex, substIndex);
        };


        // simulate   ===========================================================================
        AppUtils::SignalHandler::registerSignalHandler();
        AppUtils::Stopwatch stopWatch;
        stopWatch.start();

        reactionConditions.temperature = backgroundTemperature_K;
        for (int step=0; step<nSteps; step++) {

            double cvFieldNow_VPerM = CVFieldFct(fieldCVSetpoint_VPerM, rsSim.simulationTime());
            double svFieldNow_VPerM = SVFieldFct(fieldSVSetpoint_VPerM, rsSim.simulationTime());
            fieldMagnitude = svFieldNow_VPerM + cvFieldNow_VPerM;
            reactionConditions.electricField = fieldMagnitude;
            rsSim.performTimestep(reactionConditions, dt_s, particlesReactedFct);

            #pragma omp parallel default(none) shared(particles, mobilityScalingFct) firstprivate(mobilitScalingMode, nParticlesTotal, reducedPressure, fieldMagnitude, dt_s)
            {
                // 1 Td: = 10e-17 V*cm^2 = 10e-17 V*cm^2 = 10e-21 V*m^2
                // 2.688e19 molek/cm^3 at atmospheric pressure,
                double reducedField = (fieldMagnitude/100.0) / (reducedPressure * 2.688e2); // 2.688e19 molek/cm^3 at atmospheric pressure
                double mobilityScalingFactor = mobilityScalingFct->getInterpolatedValue(abs(reducedField));

                #pragma omp for
                for (unsigned int i = 0; i<nParticlesTotal; i++) {
                    // move particle in z axis according to particle mobility:
                    RS::ReactiveParticle* particle = particles[i].get();

                    if (mobilitScalingMode == FIELD_SCALING_FUNCTION){
                        particle->setMobility(particle->getLowFieldMobility()*mobilityScalingFactor);
                    }

                    double particleShift = particle->getMobility()*reducedPressure*fieldMagnitude*dt_s;
                    Core::Vector particleLocation = particle->getLocation();
                    particleLocation.z(particleLocation.z()+particleShift);
                    particle->setLocation(particleLocation);
                }
            }
            if (cvHighResLogPeriod > 0 && step % cvHighResLogPeriod == 0){
                cvFieldWriter->writeTimestep(cvFieldNow_VPerM, rsSim.simulationTime());
            }
            rsSim.advanceTimestep(dt_s);
            timestepWriteFct(particlesPtrs, rsSim.simulationTime(), step, false);

            //autocorrect compensation voltage, to minimize z drift (once for every single SV oscillation):
            if ( (cvMode == AUTO_CV || cvMode == MODULATED_AUTO_CV) && step % nStepsPerOscillation == 0) {
                //calculate current mean z-position:
                double buf = 0.0;
                #pragma omp parallel for reduction (+:buf) default(none) shared(particles) firstprivate(nParticlesTotal)
                for (unsigned int i = 0; i < nParticlesTotal; i++) {
                    buf += particles[i]->getLocation().z();
                }
                double currentMeanZPos = buf / nParticlesTotal;

                //update cv value:
                double diffMeanZPos = meanZPos - currentMeanZPos;
                fieldCVSetpoint_VPerM = fieldCVSetpoint_VPerM + diffMeanZPos * cvRelaxationParameter;
                cvFieldWriter->writeTimestep(std::vector<double>{fieldCVSetpoint_VPerM, currentMeanZPos}, rsSim.simulationTime());
                meanZPos = currentMeanZPos;
                logger->info("CV corrected ts:{} time:{:.2e} new CV:{} diffMeanPos:{}", step, rsSim.simulationTime(), fieldCVSetpoint_VPerM, diffMeanZPos);
            }

            if (ionsInactive>=nParticlesTotal || AppUtils::SignalHandler::isTerminationSignaled()){
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
    catch(AppUtils::TerminatedWhileCommandlineParsing& terminatedMessage){
        return terminatedMessage.returnCode();
    }
    catch(const std::invalid_argument& ia){
        std::cout << ia.what() << std::endl;
        return EXIT_FAILURE;
    }
}
