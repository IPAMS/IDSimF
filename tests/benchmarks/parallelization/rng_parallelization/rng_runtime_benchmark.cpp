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

 Simple parallel random generator runtime benchmark

 ****************************/
#include "appUtils_simulationConfiguration.hpp"
#include "appUtils_logging.hpp"
#include "appUtils_stopwatch.hpp"
#include "appUtils_signalHandler.hpp"
#include "Core_randomGenerators.hpp"
#include "CLI11.hpp"

void runParallel_oneGenerator(int nSteps, int nVals, spdlog::logger* logger){
    auto rng = Core::IndependentRandomGenerator();
    logger->info("rng address: {}", (long)&rng);

    for (int step =0; step< nSteps; ++step){
        #pragma omp parallel default(none) firstprivate(rng, nVals)
        {
            #pragma omp for
            for (int i=0; i<nVals; ++i){
                double x = rng.uniformRealRndValue();
            }
        }
    }
}

void runParallel_multipleGenerators(int nSteps, int nVals, spdlog::logger* logger){
    for (int step =0; step< nSteps; ++step){
        #pragma omp parallel default(none) firstprivate(nVals, logger)
        {
            auto rng = Core::IndependentRandomGenerator();
            //logger->info("rng address: {}", (long)&rng);
            #pragma omp for
            for (int i=0; i<nVals; ++i){
                double x = rng.uniformRealRndValue();
            }
        }
    }
}

void runParallel_randomFactory(int nSteps, int nVals, spdlog::logger* logger){
    logger->info("running runParallel_randomFactory with randomPool");
    long result = 0;
    for (int step =0; step< nSteps; ++step){
        #pragma omp parallel default(none) firstprivate(nVals, logger) shared(Core::randomGeneratorPool, result)
        {
            //logger->info("rng address: {}", (long)&rng);
            #pragma omp for
            for (int i=0; i<nVals; ++i){
                double rndVal = Core::randomGeneratorPool.getThreadElement()->normalRealRndValue();
                if (rndVal < 0.001){
                    result++;
                }
            }
        }
    }
    logger->info("result {}", result);
}

void runParallel_randomFactory_independent(int nSteps, int nVals, spdlog::logger* logger){
    logger->info("running runParallel_randomFactory_independent");
    long result = 0;
    for (int step =0; step< nSteps; ++step){
        #pragma omp parallel default(none) firstprivate(nVals, logger) shared(Core::randomGeneratorPool, result)
        {
            Core::RandomGeneratorPool::RNGPoolElement* rngElement = Core::randomGeneratorPool.getThreadElement();
            //logger->info("rng address: {}", (long)&rng);
            #pragma omp for
            for (int i=0; i<nVals; ++i){
                double rndVal = rngElement->normalRealRndValue();
                if (rndVal < 0.001){
                    result++;
                }
            }
        }
    }
    logger->info("result {}", result);
}

void runParallel_independent(int nSteps, int nVals, spdlog::logger* logger){
    logger->info("running independent");
    long result = 0;
    for (int step =0; step< nSteps; ++step){
        #pragma omp parallel default(none) firstprivate(nVals, logger) shared(Core::rdSeed, result)
        {
            std::mt19937 rngGen(Core::rdSeed());
            std::normal_distribution<double> dist;
            //logger->info("rng address: {}", (long)&rng);
            #pragma omp for
            for (int i=0; i<nVals; ++i){
                double rndVal = dist(rngGen);
                if (rndVal < 0.001){
                    result++;
                }
            }
        }
    }
    logger->info("result {}", result);
}


void runTrivial(int nSteps, int nVals, spdlog::logger* logger){
    logger->info("run trivial");
    double totalResult = 0.0;
    for (int step =0; step< nSteps; ++step){
        //logger->info("step {}", step);
        double result = 0.0;
        #pragma omp parallel default(none) firstprivate(nVals) shared(result)
        {
            #pragma omp for
            for (int i=0; i<nVals; ++i){
                double x = 0.4;
                for (int k=0; k<10000; ++k){
                    x = std::sin(x);
                }
                #pragma omp atomic
                result += x;
            }
        }
        totalResult += result;
    }
    logger->info("total result {}", totalResult);


}

int main(int argc, const char * argv[]) {

    CLI::App app{"Simple benchmark random generators", "rng benchmark"};

    int nSteps = 0;
    app.add_option("n_steps", nSteps, "number of steps")->required();

    int nVals = 0;
    app.add_option("n_vals", nVals, "number of values")->required();

    int mode = 0;
    app.add_option("mode", mode, "mode")->required();

    CLI11_PARSE(app, argc, argv);

    auto logger = AppUtils::createLogger("run.log");


    // simulate   ===========================================================================
    AppUtils::SignalHandler::registerSignalHandler();
    AppUtils::Stopwatch stopWatch;
    stopWatch.start();
    if (mode == 0){
        runParallel_multipleGenerators(nSteps, nVals, logger.get());
    }
    else if(mode== 1) {
        runParallel_oneGenerator(nSteps, nVals, logger.get());
    }
    else if(mode== 2) {
        runParallel_randomFactory(nSteps, nVals, logger.get());
    }
    else if(mode== 3) {
        runParallel_randomFactory_independent(nSteps, nVals, logger.get());
    }
    else if(mode== 4) {
        runParallel_independent(nSteps, nVals, logger.get());
    }
    else if(mode== 5){
        runTrivial(nSteps, nVals, logger.get());
    }
    stopWatch.stop();

    logger->info("----------------------");
    logger->info("CPU time: {} s", stopWatch.elapsedSecondsCPU());
    logger->info("Finished in {} seconds (wall clock time)",stopWatch.elapsedSecondsWall());

    // ======================================================================================

    return 0;
}
