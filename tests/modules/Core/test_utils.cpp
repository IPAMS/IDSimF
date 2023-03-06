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
#include <omp.h>
#include <iostream>


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

template<class GeneratorPoolType> void testGeneratorSample(
        int nSamples, expectedDistParams expectedNorm, expectedDistParams expectedUni, double margin){

    GeneratorPoolType rndGenPool;
    auto rngGen = rndGenPool.getThreadRandomSource();

    std::vector<double> normalVals = std::vector<double>();
    std::vector<double> uniformVals = std::vector<double>();

    for (int i=0; i<nSamples;i++){
        normalVals.push_back(rngGen->normalRealRndValue());
        uniformVals.push_back(rngGen->uniformRealRndValue());
    }

    randomSampleParams srNorm = calculateParamsFromSample(normalVals);
    randomSampleParams srUni = calculateParamsFromSample(uniformVals);


    checkDistParams(srNorm,expectedNorm,margin);
    checkDistParams(srUni,expectedUni,margin);
}

template<class GeneratorPoolType> void testUniformCustomDistribution(int nSamples, double min, double max){
    GeneratorPoolType rndGenPool;
    Core::RndDistPtr uniDist = rndGenPool.getUniformDistribution(min, max);
    std::vector<double> uniformVals = std::vector<double>();

    for (int i=0; i<nSamples; ++i){
        uniformVals.push_back(uniDist->rndValue());
        REQUIRE(uniformVals.back() > min);
        REQUIRE(uniformVals.back() < max);
    }
}

TEST_CASE("Test random bit sources") {
    SECTION("Test Mersenne Bit Source"){
        std::vector<int> testVector1 = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
        std::vector<int> testVector2 = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};

        SECTION("Test with truly random seed"){

            Core::MersenneBitSource mersenneSource1;
            Core::MersenneBitSource mersenneSource2;

            std::shuffle(testVector1.begin(), testVector1.end(), mersenneSource1);
            std::shuffle(testVector2.begin(), testVector2.end(), mersenneSource2);
            CHECK( (testVector1[0] != 1 && testVector1[1] != 2) );
            CHECK( testVector1[0] != testVector2[0]);
            CHECK( testVector1[1] != testVector2[1]);
        }

        SECTION("Test reproductivity with constant seed"){

            Core::MersenneBitSource mersenneSource1;
            Core::MersenneBitSource mersenneSource2;

            mersenneSource1.seed(200);
            mersenneSource2.seed(200);

            std::shuffle(testVector1.begin(), testVector1.end(), mersenneSource1);
            std::shuffle(testVector2.begin(), testVector2.end(), mersenneSource2);
            CHECK( (testVector1[0] != 1 && testVector1[1] != 2) );
            CHECK( testVector1[0] == testVector2[0]);
            CHECK( testVector1[1] == testVector2[1]);
            CHECK( testVector1[3] == testVector2[3]);

        }

        SECTION("Test different constant seed"){

            Core::MersenneBitSource mersenneSource1;
            Core::MersenneBitSource mersenneSource2;

            mersenneSource1.seed(200);
            mersenneSource2.seed(300);

            std::shuffle(testVector1.begin(), testVector1.end(), mersenneSource1);
            std::shuffle(testVector2.begin(), testVector2.end(), mersenneSource2);
            CHECK( (testVector1[0] != 1 && testVector1[1] != 2) );
            CHECK( testVector1[0] != testVector2[0]);
            CHECK( testVector1[1] != testVector2[1]);
            CHECK( testVector1[3] != testVector2[3]);

        }

    }


    SECTION("Test Test Bit Source"){
        Core::TestBitSource testSource1;
        Core::TestBitSource testSource2;

        std::vector<int> testVector1 = {1,2,3,4,5,6,7,8};
        std::vector<int> testVector2 = {1,2,3,4,5,6,7,8};

        std::shuffle(testVector1.begin(), testVector1.end(), testSource1);
        std::shuffle(testVector2.begin(), testVector2.end(), testSource2);

        CHECK( (testVector1[0] != 1 && testVector1[1] != 2) );
        CHECK(testVector1 == testVector2);
        CHECK(testSource1() != testSource1());
    }

    SECTION("Test SplitMix64 reproductivity with constant seed"){
            std::vector<int> testVector1 = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
            std::vector<int> testVector2 = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};

            Core::SplitMix64TestBitSource splitMix64Source1 = Core::SplitMix64TestBitSource();
            Core::SplitMix64TestBitSource splitMix64Source2 = Core::SplitMix64TestBitSource();

            std::shuffle(testVector1.begin(), testVector1.end(), splitMix64Source1);
            std::shuffle(testVector2.begin(), testVector2.end(), splitMix64Source2);
            CHECK( (testVector1[0] != 1 && testVector1[1] != 2) );
            CHECK( testVector1[0] == testVector2[0]);
            CHECK( testVector1[1] == testVector2[1]);
            CHECK( testVector1[3] == testVector2[3]);
    }

    SECTION("Test Xoshiro reproductivity with constant seed"){
            std::vector<int> testVector1 = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
            std::vector<int> testVector2 = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};

            Core::Xoshiro256pTestBitSource xoshiro256pSource1 = Core::Xoshiro256pTestBitSource();
            Core::Xoshiro256pTestBitSource xoshiro256pSource2 = Core::Xoshiro256pTestBitSource();

            std::shuffle(testVector1.begin(), testVector1.end(), xoshiro256pSource1);
            std::shuffle(testVector2.begin(), testVector2.end(), xoshiro256pSource2);
            CHECK( (testVector1[0] != 1 && testVector1[1] != 2) );
            CHECK( testVector1[0] == testVector2[0]);
            CHECK( testVector1[1] == testVector2[1]);
            CHECK( testVector1[3] == testVector2[3]);
    }

}

TEST_CASE("Test productive random generator pool") {

    int nMaxThreads = omp_get_max_threads();
    Core::RandomGeneratorPool rngPool;

    SECTION("Generators in rng generator pool should produce independent randomness"){
        if (nMaxThreads > 1){
            auto rng0 = rngPool.getRandomSource(0);
            auto rng1 = rngPool.getRandomSource(1);

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
            CHECK(rng0->uniformRealRndValue() != rng1->uniformRealRndValue());
#pragma GCC diagnostic pop

        }
    }

    SECTION("There should be a rng for all threads"){
        CHECK_NOTHROW(Core::globalRandomGeneratorPool->getRandomSource(static_cast<std::size_t>(nMaxThreads-1)));
    }
}


TEST_CASE( "Test productive random distributions", "[Core][random]") {

    SECTION("Uniform random distribution should have the correct mean and deviation") {
        Core::RandomGeneratorPool rngPool;
        auto mersenneBitSource = rngPool.getThreadRandomSource()->getRandomBitSource();

        int nSamples = 10000;

        Core::UniformRandomDistribution dist(0.0, 1.0, mersenneBitSource);

        std::vector<double> vals;
        for (int i=0; i<nSamples;i++){
            vals.push_back(dist.rndValue());
        }
        randomSampleParams sr = calculateParamsFromSample(vals);

        bool meanInRage = ((sr.mean<0.51) && (sr.mean>0.49));
        CHECK(meanInRage);
        bool stdDevInRage = ((sr.stdDev<0.30) && (sr.stdDev>0.26));
        CHECK(stdDevInRage);
    }

    SECTION("Productive rng should produce correct uniform and normal distributions"){
        int nSamples = 10000;
        expectedDistParams expectedNorm{0.0, 1.0};
        expectedDistParams expectedUni{0.5, 0.28};

        testGeneratorSample<Core::RandomGeneratorPool>(nSamples, expectedNorm, expectedUni, 0.02);
        testUniformCustomDistribution<Core::RandomGeneratorPool>(1000, 2.0, 6.0);
    }

    SECTION("Constant seeding of productive random generator pool elements should be possible") {
        Core::RandomGeneratorPool rngPool1;
        Core::RandomGeneratorPool rngPool2;

        int nSamples = 10;

        for (int i=0; i<nSamples; ++i){
            CHECK(rngPool1.getThreadRandomSource()->uniformRealRndValue() != Approx(rngPool2.getThreadRandomSource()->uniformRealRndValue()));
        }

        rngPool1.setSeedForElements(100);
        rngPool2.setSeedForElements(100);

        for (int i=0; i<nSamples; ++i){
            CHECK(rngPool1.getThreadRandomSource()->uniformRealRndValue() == Approx(rngPool2.getThreadRandomSource()->uniformRealRndValue()));
        }

    }
}

TEST_CASE( "Test testing random distributions", "[Core][random]") {

    Core::TestRandomGeneratorPool rngPool;
    Core::TestRandomGeneratorPool::TestRNGPoolElement* rngPoolElem = rngPool.getThreadRandomSource();

    std::vector<double> vals;

    SECTION("Test random generator should generate predefined deterministic normal random samples") {

        for (int i = 0; i<10; ++i) {
            vals.push_back(rngPoolElem->normalRealRndValue());
        }

        REQUIRE(Approx(vals[2])==0.046029231875);
        REQUIRE(Approx(vals[4])==0.176058757790);
        REQUIRE(Approx(vals[8])==-0.628735061607);
    }

    SECTION("Test random generator should generate predefined deterministic uniform random samples") {

        for (int i = 0; i<10; ++i) {
            vals.push_back(rngPoolElem->uniformRealRndValue());
        }

        REQUIRE(Approx(vals[1])==0.900613166242);
        REQUIRE(Approx(vals[3])==0.498839039733);
        REQUIRE(Approx(vals[7])==0.437690201598);
    }

    SECTION("Test random generator produces skewed distribution"){

        int nSamples = 1000;
        expectedDistParams expectedNorm{0.1615938885, 0.8400041942};
        expectedDistParams expectedUni{0.5253988, 0.2888807331};

        testGeneratorSample<Core::TestRandomGeneratorPool>(nSamples, expectedNorm, expectedUni, 1e-7);

        testUniformCustomDistribution<Core::TestRandomGeneratorPool>(1000, 2.0, 6.0);
    }
}

TEST_CASE( "Test xoshiro256+ test random distributions (constant seeding)", "[Core][random]") {

    Core::XoshiroTestRandomGeneratorPool rngPool;
    Core::XoshiroTestRandomGeneratorPool::XoshiroTestRNGPoolElement* rngPoolElem = rngPool.getThreadRandomSource();

    std::vector<double> vals;

    SECTION("Test random generator should generate predefined deterministic normal random samples") {

        for (int i = 0; i<10; ++i) {
            vals.push_back(rngPoolElem->normalRealRndValue());
        }

        REQUIRE(Approx(vals[2])==-1.4300851003);
        REQUIRE(Approx(vals[4])==-0.8904677016);
        REQUIRE(Approx(vals[8])==-2.7963889313);
    }

    SECTION("Test random generator should generate predefined deterministic uniform random samples") {

        for (int i = 0; i<10; ++i) {
            vals.push_back(rngPoolElem->uniformRealRndValue());
        }

        REQUIRE(Approx(vals[1])==0.9264488555);
        REQUIRE(Approx(vals[3])==0.2340490438);
        REQUIRE(Approx(vals[7])==0.9594879628);
    }

     SECTION("Uniform random distribution should have the correct mean and deviation") {
        auto xoshiroBitSource = rngPool.getThreadRandomSource()->getRandomBitSource();

        int nSamples = 10000;

        Core::UniformTestDistributionXoshiro dist(0.0, 1.0, xoshiroBitSource);

        std::vector<double> vals;
        for (int i=0; i<nSamples;i++){
            vals.push_back(dist.rndValue());
        }
        randomSampleParams sr = calculateParamsFromSample(vals);

        bool meanInRage = ((sr.mean<0.51) && (sr.mean>0.49));
        CHECK(meanInRage);
        bool stdDevInRage = ((sr.stdDev<0.30) && (sr.stdDev>0.26));
        CHECK(stdDevInRage);
    }

    SECTION("Test rng should produce correct uniform and normal distributions"){
        int nSamples = 10000;
        expectedDistParams expectedNorm{0.0, 1.0};
        expectedDistParams expectedUni{0.5, 0.28};

        testGeneratorSample<Core::XoshiroTestRandomGeneratorPool>(nSamples, expectedNorm, expectedUni, 0.02);
        testUniformCustomDistribution<Core::XoshiroTestRandomGeneratorPool>(1000, 2.0, 6.0);
    }
}
