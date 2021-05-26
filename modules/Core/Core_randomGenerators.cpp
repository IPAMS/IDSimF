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
 ****************************/

#include "Core_randomGenerators.hpp"

std::random_device Core::rdSeed; //seed generator
std::mt19937 Core::internalRNG(Core::rdSeed()); //internal real random generator

std::unique_ptr<Core::AbstractRandomGenerator> Core::globalRandomGenerator =
        std::make_unique<Core::RandomGenerator>();

/**
 * Constructs a uniform random distribution in the interval [0.0, 1.0]
 */
Core::UniformRandomDistribution::UniformRandomDistribution(): internalUniformDist_(0.0, 1.0)
{}

/**
 * Construct a custom uniform random distribution with a custom interval
 * @param min lower boundary of the interval
 * @param max upper boundary of the interval
 */
Core::UniformRandomDistribution::UniformRandomDistribution(double min, double max): internalUniformDist_(min, max)
{}

/**
 * Generate random value in the uniform distribution
 * @return uniformly distributed random value
 */
double Core::UniformRandomDistribution::rndValue() {
    return internalUniformDist_(Core::internalRNG);
}

/**
 * Construct normal random distribution
 */
Core::NormalRandomDistribution::NormalRandomDistribution():
    internalNormalDist_()
{}

/**
 * Generate random value
 * @return normal distributed random value
 */
double Core::NormalRandomDistribution::rndValue() {
    return internalNormalDist_(Core::internalRNG);
}

/**
 * Constructs a test distribution in the interval [0.0, 1.0]
 */
Core::UniformTestDistribution::UniformTestDistribution():
    sampleIndex_(0),
    min_(0),
    interval_(1.0)
{}

/**
 * Construct a custom test distribution with a custom interval
 * @param min lower boundary of the interval
 * @param max upper boundary of the interval
 */
Core::UniformTestDistribution::UniformTestDistribution(double min, double max):
    sampleIndex_(0),
    min_(min),
    interval_(max-min)
{}

/**
 * Generate non random test value in the uniform distribution
 * @return uniformly distributed test value
 */
double Core::UniformTestDistribution::rndValue() {
    sampleIndex_ = (sampleIndex_+1) % Core::UNIFORM_TEST_SAMPLES.size();
    return min_ + Core::UNIFORM_TEST_SAMPLES[sampleIndex_]*interval_;
}

/**
 * Construct normal test distribution
 */
Core::NormalTestDistribution::NormalTestDistribution(): sampleIndex_(0)
{}

/**
 * Generate normal test value
 * @return non random normal distributed test value
 */
double Core::NormalTestDistribution::rndValue() {
    sampleIndex_ = (sampleIndex_+1) % Core::NORMAL_TEST_SAMPLES.size();
    return Core::NORMAL_TEST_SAMPLES[sampleIndex_];
}

/**
 * Construct production use RandomGenerator
 */
Core::RandomGenerator::RandomGenerator():
    uniformDistribution_(),
    normalDistribution_()
{}

/**
 * Generate uniformly distributed in the interval [0.0, 1.0]
 * @return uniformly distributed random value in [0.0, 1.0]
 */
double Core::RandomGenerator::uniformRealRndValue() {
    return uniformDistribution_.rndValue();
}

/**
 * Generate normal distributed value
 * @return normal distributed random value
 */
double Core::RandomGenerator::normalRealRndValue() {
    return normalDistribution_.rndValue();
}

/**
 * Generate a custom uniform random distribution in the interval [min, max]
 * @param min lower boundary of the random distribution
 * @param max upper boundary of the random distribution
 * @return uniform random distribution in the interval [min, max]
 */
std::unique_ptr<Core::RandomDistribution> Core::RandomGenerator::getUniformDistribution(double min, double max) {
    return std::make_unique<Core::UniformRandomDistribution>(min, max);
}


/**
 * Construct test value random generator
 */
Core::TestRandomGenerator::TestRandomGenerator():
    uniformDistribution_(),
    normalDistribution_()
{}

/**
 * Generate uniformly distributed *non* random test value in the interval [0.0, 1.0]
 * @return uniformly distributed test value
 */
double Core::TestRandomGenerator::uniformRealRndValue() {
    return uniformDistribution_.rndValue();
}

/**
 * Generate normal distributed *non* random test value
 * @return normal distributed test value
 */
double Core::TestRandomGenerator::normalRealRndValue() {
    return normalDistribution_.rndValue();
}

/**
 * Generate a custom uniform *non* random test distribution in the interval [min, max]
 * @param min lower boundary of the random distribution
 * @param max upper boundary of the random distribution
 * @return uniform test distribution in the interval [min, max]
 */
std::unique_ptr<Core::RandomDistribution> Core::TestRandomGenerator::getUniformDistribution(double min, double max) {
    return std::make_unique<Core::UniformTestDistribution>(min, max);
}
