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
#include <omp.h>

std::random_device Core::rdSeed; //seed generator

std::unique_ptr<Core::AbstractRandomGenerator> Core::globalRandomGenerator =
        std::make_unique<Core::PooledRandomGenerator>();

std::unique_ptr<Core::AbstractRandomGeneratorPool> Core::globalRandomGeneratorPool =
        std::make_unique<Core::RandomGeneratorPool>();


Core::MersenneBitSource::MersenneBitSource():
randomSource(Core::rdSeed())
{}

std::mt19937::result_type Core::MersenneBitSource::max() {
    return randomSource.max();
}

std::mt19937::result_type Core::MersenneBitSource::min() {
    return randomSource.min();
}

std::mt19937::result_type Core::MersenneBitSource::operator()() {
    return randomSource();
}

Core::TestBitSource::TestBitSource():
sampleIndex_(0)
{}

unsigned int Core::TestBitSource::min() {
    return 0;
}

unsigned int Core::TestBitSource::max() {
    return 4294967295;
}

unsigned int Core::TestBitSource::operator()() {
    sampleIndex_ = (sampleIndex_+1) % Core::UNIFORM_RANDOM_BITS.size();
    return Core::UNIFORM_RANDOM_BITS[sampleIndex_];
}

Core::RandomGeneratorPool::RNGPoolElement::RNGPoolElement():
rngGenerator_(Core::rdSeed())
{}

double Core::RandomGeneratorPool::RNGPoolElement::uniformRealRndValue() {
    return uniformDist_(rngGenerator_);
}

double Core::RandomGeneratorPool::RNGPoolElement::normalRealRndValue() {
    return normalDist_(rngGenerator_);
}

std::mt19937* Core::RandomGeneratorPool::RNGPoolElement::getRNG() {
    return &rngGenerator_;
}

/**
 * Constructs the internal random number generator
 */
Core::RandomGeneratorPool::RandomGeneratorPool() {
    int nMaxThreads_ = omp_get_max_threads();
    for (int i=0; i<nMaxThreads_; ++i){
        elements_.emplace_back(std::make_unique<Core::RandomGeneratorPool::RNGPoolElement>());
    }
}

Core::RandomGeneratorPool::RNGPoolElement* Core::RandomGeneratorPool::getThreadElement() {
    return elements_[static_cast<std::size_t>(omp_get_thread_num())].get();
}

Core::RandomGeneratorPool::RNGPoolElement* Core::RandomGeneratorPool::getElement(std::size_t index) {
    return elements_.at(index).get();
}


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
    return 0.0; //internalUniformDist_(Core::globalRandomGeneratorPool->getThreadElement().);
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
    return 0.0;// internalNormalDist_(Core::oldInternalRNG);
}

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
Core::PooledRandomGenerator::PooledRandomGenerator()
{}

/**
 * Generate uniformly distributed in the interval [0.0, 1.0]
 * @return uniformly distributed random value in [0.0, 1.0]
 */
double Core::PooledRandomGenerator::uniformRealRndValue() {
    return Core::globalRandomGeneratorPool->getThreadElement()->uniformRealRndValue();
}

/**
 * Generate normal distributed value
 * @return normal distributed random value
 */
double Core::PooledRandomGenerator::normalRealRndValue() {
    return Core::globalRandomGeneratorPool->getThreadElement()->normalRealRndValue();
}

/**
 * Generate a custom uniform random distribution in the interval [min, max]
 * @param min lower boundary of the random distribution
 * @param max upper boundary of the random distribution
 * @return uniform random distribution in the interval [min, max]
 */
std::unique_ptr<Core::RandomDistribution> Core::PooledRandomGenerator::getUniformDistribution(double min, double max) {
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
