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

std::unique_ptr<Core::AbstractRandomGeneratorPool> Core::globalRandomGeneratorPool =
        std::make_unique<Core::RandomGeneratorPool>();


Core::MersenneBitSource::MersenneBitSource():
internalRandomSource(Core::rdSeed())
{}

Core::rndBit_type Core::MersenneBitSource::max() {
    return internalRandomSource.max();
}

Core::rndBit_type Core::MersenneBitSource::min() {
    return internalRandomSource.min();
}

Core::rndBit_type Core::MersenneBitSource::operator()() {
    return internalRandomSource();
}

Core::TestBitSource::TestBitSource():
sampleIndex_(0)
{}

Core::rndBit_type Core::TestBitSource::min() {
    return 0;
}

Core::rndBit_type Core::TestBitSource::max() {
    return 4294967295;
}

Core::rndBit_type Core::TestBitSource::operator()() {
    sampleIndex_ = (sampleIndex_+1) % Core::UNIFORM_RANDOM_BITS.size();
    return Core::UNIFORM_RANDOM_BITS[sampleIndex_];
}


/**
 * Construct a custom uniform random distribution with a custom interval
 * @param min lower boundary of the interval
 * @param max upper boundary of the interval
 */
Core::UniformRandomDistribution::UniformRandomDistribution(double min, double max, Core::RandomBitSource<rndBit_type>* randomSource):
randomSource_(randomSource),
internalUniformDist_(min, max)
{}

/**
 * Generate random value in the uniform distribution
 * @return uniformly distributed random value
 */
double Core::UniformRandomDistribution::rndValue(){
    return internalUniformDist_(*randomSource_);
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


double Core::RandomGeneratorPool::RNGPoolElement::uniformRealRndValue() {
    return uniformDist_(rngGenerator_.internalRandomSource);
}

double Core::RandomGeneratorPool::RNGPoolElement::normalRealRndValue() {
    return normalDist_(rngGenerator_.internalRandomSource);
}

Core::MersenneBitSource* Core::RandomGeneratorPool::RNGPoolElement::getRandomBitSource() {
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

Core::RndDistPtr Core::RandomGeneratorPool::getUniformDistribution(double min, double max) {
    return std::make_unique<Core::UniformRandomDistribution>(min, max, this->getThreadElement()->getRandomBitSource());
}

Core::RandomGeneratorPool::RNGPoolElement* Core::RandomGeneratorPool::getThreadElement() {
    return elements_[static_cast<std::size_t>(omp_get_thread_num())].get();
}

Core::RandomGeneratorPool::RNGPoolElement* Core::RandomGeneratorPool::getElement(std::size_t index) {
    return elements_.at(index).get();
}

double Core::TestRandomGeneratorPool::TestPoolElement::uniformRealRndValue() {
    return uniformDist_.rndValue();
}

double Core::TestRandomGeneratorPool::TestPoolElement::normalRealRndValue() {
    return normalDist_.rndValue();
}

Core::TestBitSource* Core::TestRandomGeneratorPool::TestPoolElement::getRandomBitSource() {
    return &rngGenerator_;
}

Core::RndDistPtr Core::TestRandomGeneratorPool::getUniformDistribution(double min, double max) {
    return std::make_unique<Core::UniformTestDistribution>(min, max);
}

Core::TestRandomGeneratorPool::TestPoolElement* Core::TestRandomGeneratorPool::getThreadElement() {
    return &element_;
}

Core::TestRandomGeneratorPool::TestPoolElement* Core::TestRandomGeneratorPool::getElement(std::size_t) {
    return &element_;
}
