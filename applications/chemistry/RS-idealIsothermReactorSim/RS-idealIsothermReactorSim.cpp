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
 RS-idealIsothermReactorSim.cpp

 Simple chemical kinetics in an isotherm ideally mixed reactor

 ****************************/
#include "Core_randomGenerators.hpp"
#include "RS_Simulation.hpp"
#include "RS_SimulationConfiguration.hpp"
#include "RS_ConfigFileParser.hpp"
#include "RS_ConcentrationFileWriter.hpp"
#include "appUtils_simulationConfiguration.hpp"
#include "appUtils_logging.hpp"
#include "appUtils_stopwatch.hpp"
#include "appUtils_signalHandler.hpp"
#include "appUtils_commandlineParser.hpp"
#include <iostream>
#include <cmath>

int main(int argc, const char * argv[]) {

    try{
        // parse commandline / create conf and logger ===================================================
        AppUtils::CommandlineParser cmdLineParser(argc, argv, "BT-idealIsothermReactorSim",
                "RS Simulation in an ideally mixed, isotherm reactor", true);
        std::string resultFilename = cmdLineParser.projectName() + "_result.txt";
        auto logger = cmdLineParser.logger();
        auto simConf = cmdLineParser.simulationConfiguration();

        std::string rsConfigFileName = simConf->pathRelativeToConfFile(
                simConf->stringParameter("reaction_configuration"));
        RS::ConfigFileParser parser = RS::ConfigFileParser();
        RS::Simulation sim = RS::Simulation(parser.parseFile(rsConfigFileName));
        RS::SimulationConfiguration* rsSimConf = sim.simulationConfiguration();

        std::vector<int> nParticles = simConf->intVectorParameter("n_particles");
        int nSteps = simConf->intParameter("sim_time_steps");
        int concentrationWriteInterval = simConf->intParameter("concentrations_write_interval");
        double dt_s = simConf->doubleParameter("dt_s");
        double backgroundTemperature_K = simConf->doubleParameter("background_temperature_K");

        RS::ConcentrationFileWriter resultFilewriter(resultFilename);
        // ======================================================================================

        // init simulation  =====================================================================

        // create and add simulation particles:
        std::size_t nParticlesTotal = 0;
        std::vector<uniqueReactivePartPtr> particles;
        for (std::size_t i=0; i<nParticles.size(); ++i) {
            RS::Substance *subst = rsSimConf->getAllDiscreteSubstances().at(i);
            for (int k = 0; k < nParticles[i]; ++k) {
                uniqueReactivePartPtr particle = std::make_unique<RS::ReactiveParticle>(subst);
                sim.addParticle(particle.get(), nParticlesTotal);
                particles.push_back(std::move(particle));
                nParticlesTotal++;
            }
        }
        RS::ReactionConditions reactionConditions = RS::ReactionConditions();
        reactionConditions.temperature = backgroundTemperature_K;
        reactionConditions.electricField = 0.0;
        reactionConditions.pressure = 0.0;

        resultFilewriter.initFile(rsSimConf);
        // ======================================================================================


        // simulate   ===========================================================================
        AppUtils::SignalHandler::registerSignalHandler();
        AppUtils::Stopwatch stopWatch;
        stopWatch.start();
        for (int step=0; step<nSteps; ++step) {
            if (AppUtils::SignalHandler::isTerminationSignaled()){ // terminate simulation loop if termination was signaled
                break;
            }
            if (step % concentrationWriteInterval ==0) {
                sim.logConcentrations(logger);
                resultFilewriter.writeTimestep(sim);
            }

            sim.performTimestep(reactionConditions, dt_s);
            sim.advanceTimestep(dt_s);
        }
        resultFilewriter.closeFile();
        stopWatch.stop();

        logger->info("----------------------");
        logger->info("Reaction Events:");
        sim.logReactionStatistics(logger);
        logger->info("----------------------");
        logger->info("total reaction events: {} ill events: {}", sim.totalReactionEvents(), sim.illEvents());
        logger->info("ill fraction: {}", sim.illEvents() / (double) sim.totalReactionEvents());

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
