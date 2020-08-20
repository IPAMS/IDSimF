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
 test_utils.cpp

 Testing of utilities from core module (primarily random number generators)

 ****************************/

#include "Core_constants.hpp"
#include "Core_randomGenerators.hpp"
#include "catch.hpp"


struct randomSampleParams{
    double sum;
    double mean;
    double stdDev;
};

struct expectedDistParams{
    double mean;
    double stdDev;
};

randomSampleParams calculateParamsFromSample(std::vector<double> sample){
    double sum = 0.0;
    for (double smp : sample){
        sum = sum + smp;
    }

    double mean = sum / sample.size();
    double diffSum = 0.0;
    for (size_t i=0; i < sample.size(); ++i){
        diffSum = diffSum + std::pow(sample[i]-mean,2.0);
    }
    double stdDev = std::sqrt(1.0 / (sample.size()-1.0) * diffSum);

    randomSampleParams result{sum, mean, stdDev};

    return result;
}

template<class DistType> randomSampleParams generateDistributionSample(int nSamples){
    DistType dist;
    std::vector<double> vals = std::vector<double>();
    for (int i=0; i<nSamples;i++){
        vals.push_back(dist.rndValue());
    }

    randomSampleParams result = calculateParamsFromSample(vals);

    return result;
}

void checkDistParams(randomSampleParams rsParam, expectedDistParams expectedParam, double margin){
    REQUIRE(rsParam.mean == Approx(expectedParam.mean).margin(margin));
    REQUIRE(rsParam.stdDev == Approx(expectedParam.stdDev).margin(margin));
}

template<class GeneratorType> void testGeneratorSample(
        int nSamples, expectedDistParams expectedNorm, expectedDistParams expectedUni, double margin){

    GeneratorType rndGen;

    std::vector<double> normalVals = std::vector<double>();
    std::vector<double> uniformVals = std::vector<double>();

    for (int i=0; i<nSamples;i++){
        normalVals.push_back(rndGen.normalRealRndValue());
        uniformVals.push_back(rndGen.uniformRealRndValue());
    }

    randomSampleParams srNorm = calculateParamsFromSample(normalVals);
    randomSampleParams srUni = calculateParamsFromSample(uniformVals);


    checkDistParams(srNorm,expectedNorm,margin);
    checkDistParams(srUni,expectedUni,margin);
}

template<class GeneratorType> void testUniformCustomDistribution(int nSamples, double min, double max){
    GeneratorType rndGen;
    Core::RndDistPtr uniDist = rndGen.getUniformDistribution(min, max);
    std::vector<double> uniformVals = std::vector<double>();

    for (int i=0; i<nSamples;i++){
        uniformVals.push_back(uniDist->rndValue());
        REQUIRE(uniformVals.back() > min);
        REQUIRE(uniformVals.back() < max);
    }
}


TEST_CASE( "Test random distributions", "[Core][random]") {

    SECTION("Uniform random distribution should have the correct mean and deviation") {

        int nSamples = 10000;
        randomSampleParams sr = generateDistributionSample<Core::UniformRandomDistribution>(nSamples);

        bool meanInRage = ((sr.mean<0.51) && (sr.mean>0.49));
        REQUIRE(meanInRage);
        bool stdDevInRage = ((sr.stdDev<0.30) && (sr.stdDev>0.26));
        REQUIRE(stdDevInRage);
    }

    SECTION("Normal random distribution should have the correct mean and deviation") {
        int nSamples = 1000;
        randomSampleParams sr = generateDistributionSample<Core::NormalRandomDistribution>(nSamples);

        REQUIRE(sr.mean<0.1);
        bool stdDevInRage = ((sr.stdDev<1.05) && (sr.stdDev>0.95));
        REQUIRE(stdDevInRage);
    }
}

TEST_CASE( "Test random generators", "[Core][random]") {
    SECTION("Production random generator should generate correct uniform random samples") {
        int nSamples = 10000;
        expectedDistParams expectedNorm{0.0, 1.0};
        expectedDistParams expectedUni{0.5, 0.28};

        testGeneratorSample<Core::RandomGenerator>(nSamples, expectedNorm, expectedUni, 0.02);

        testUniformCustomDistribution<Core::RandomGenerator>(1000, 2.0, 6.0);
    }

    SECTION("Test deterministic test random generators"){

        SECTION("Test individual test generator samples") {

            Core::TestRandomGenerator testRng;
            std::vector<double> vals;

            SECTION("Test random generator should generate predefined deterministic normal random samples") {

                for (int i = 0; i<10; ++i) {
                    vals.push_back(testRng.normalRealRndValue());
                }

                REQUIRE(Approx(vals[2])==0.046029231875);
                REQUIRE(Approx(vals[4])==0.176058757790);
                REQUIRE(Approx(vals[8])==-0.628735061607);
            }

            SECTION("Test random generator should generate predefined deterministic uniform random samples") {

                for (int i = 0; i<10; ++i) {
                    vals.push_back(testRng.uniformRealRndValue());
                }

                REQUIRE(Approx(vals[1])==0.900613166242);
                REQUIRE(Approx(vals[3])==0.498839039733);
                REQUIRE(Approx(vals[7])==0.437690201598);
            }
        }

        SECTION("Test random generator produces skewed disitribution"){

            int nSamples = 1000;
            expectedDistParams expectedNorm{0.1615938885, 0.8400041942};
            expectedDistParams expectedUni{0.5253988, 0.2888807331};

            testGeneratorSample<Core::TestRandomGenerator>(nSamples, expectedNorm, expectedUni, 1e-7);

            testUniformCustomDistribution<Core::TestRandomGenerator>(1000, 2.0, 6.0);
        }

    }

}