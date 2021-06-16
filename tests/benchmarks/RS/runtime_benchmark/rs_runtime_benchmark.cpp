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

 A simple runtime benchmarks for RS for benchmarking and profiling

 ***************************/

#include "RS_Simulation.hpp"
#include "RS_SimulationConfiguration.hpp"
#include "RS_ConfigFileParser.hpp"
#include "RS_ConcentrationFileWriter.hpp"
#include "appUtils_stopwatch.hpp"
#include <iostream>
#include <cmath>

int main() {

    //read and prepare chemical configuration ===============================================
    std::string confFileName = "RS_waterCl_static.conf";
    RS::ConfigFileParser parser = RS::ConfigFileParser();
    RS::Simulation rsSim = RS::Simulation(parser.parseFile(confFileName));
    RS::SimulationConfiguration* simConf = rsSim.simulationConfiguration();

    // init simulation  =====================================================================

    // create and add simulation particles:
    int nParticles = 10000;
    std::vector<uniqueReactivePartPtr>particles;
    std::vector<BTree::Particle*>particlesPtrs;
    //std::vector<std::vector<double>> trajectoryAdditionalParams;

    RS::Substance *subst = simConf->substance(0);
    for (int i=0; i<nParticles;i++) {

        uniqueReactivePartPtr particle = std::make_unique<RS::ReactiveParticle>(subst);
        particle->setLocation({0.01*i, 0.0, 0.0});
        particlesPtrs.push_back(particle.get());
        rsSim.addParticle(particle.get(), i);
        particles.push_back(std::move(particle));

    }

    RS::ReactionConditions reactionConditions = RS::ReactionConditions();
    reactionConditions.temperature = 300;
    reactionConditions.pressure = 100000;
    reactionConditions.electricField = 0.0;
    // ======================================================================================


    // simulate   ===========================================================================
    int nSteps = 10000;
    int concentrationPrintInterval = 100;
    double dt_s = 5e-11;

    AppUtils::Stopwatch stopWatch;
    stopWatch.start();

    for (int step=0; step<nSteps; step++) {

        for (int i = 0; i < nParticles; i++) {
            rsSim.react(i, reactionConditions, dt_s);
        }
        rsSim.advanceTimestep(dt_s);

        if (step % concentrationPrintInterval ==0) {
            rsSim.printConcentrations();
        }

    }

    stopWatch.stop();

    std::cout << "reaction events:" << rsSim.totalReactionEvents() << " ill events:" << rsSim.illEvents() << std::endl;
    std::cout << "ill fraction: " << rsSim.illEvents() / (double) rsSim.totalReactionEvents() << std::endl;
    std::cout << "elapsed wall time:"<< stopWatch.elapsedSecondsWall()<<std::endl;
    std::cout << "elapsed cpu time:"<< stopWatch.elapsedSecondsCPU()<<std::endl;

    // ======================================================================================

    return 0;
}
