#include "CollisionModel_MathFunctions.hpp"
#include "appUtils_stopwatch.hpp"
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>

double generateFct()
{
    static double v = 1;
    v = v + 0.01;
    return v;
}

int main() {
    std::cout << "Benchmark errorfunction" << std::endl;
    AppUtils::Stopwatch stopWatch;
    stopWatch.start();

    unsigned int nSamples = 80000000;
    std::vector<double> xVec(nSamples);
    std::generate(xVec.begin(), xVec.end(), generateFct);
    std::vector<double> result_std(nSamples);

    for (unsigned int i=0; i < nSamples; ++i){
        result_std[i] = std::erf(xVec[i]);
    }

    stopWatch.stop();

    std::cout << "elapsed wall time:"<< stopWatch.elapsedSecondsWall()<<std::endl;
    std::cout << "elapsed cpu time:"<< stopWatch.elapsedSecondsCPU()<<std::endl;
    std::cout << "result[0]: " << result_std[0] << std::endl;
    return 0;
}
