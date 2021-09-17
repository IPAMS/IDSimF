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

 Core_randomGenerators.hpp

 Random number generators for the Core / IDSimF Project

 Simple abstraction classes for random number generation,
 for simulation and testing use

 ****************************/

#ifndef BTree_randomGenerators
#define BTree_randomGenerators

#include "Core_randomTestSamples.hpp"
#include <random>
#include <vector>
#include <memory>

namespace Core{

    using rndBit_type = std::mt19937::result_type;

    extern std::random_device rdSeed; ///< global seed generator

    template <class result_T>
    class RandomBitSource{
    public:
        typedef result_T result_type;
        virtual ~RandomBitSource() =default;
        virtual result_T min() =0;
        virtual result_T max() =0;
        virtual result_T operator()() =0;
    };

    class MersenneBitSource: public RandomBitSource<rndBit_type>{
    public:
        MersenneBitSource();
        rndBit_type min() override;
        rndBit_type max() override;
        rndBit_type operator()() override;

        std::mt19937 internalRandomSource;
    };

    class TestBitSource: public RandomBitSource<rndBit_type>{
    public:
        TestBitSource();
        rndBit_type min() override;
        rndBit_type max() override;
        rndBit_type operator()() override;

    private:
        std::size_t sampleIndex_;
    };

    class AbstractRNGPoolElement{
    public:
        virtual ~AbstractRNGPoolElement() =default;
        virtual double uniformRealRndValue() =0;
        virtual double normalRealRndValue() =0;
        virtual RandomBitSource<rndBit_type>* getRandomBitSource() =0;
    };

    class AbstractRandomGeneratorPool{
    public:
        virtual ~AbstractRandomGeneratorPool() =default;
//        virtual RndDistPtr getUniformDistribution(double min, double max) =0;
        virtual AbstractRNGPoolElement* getThreadElement() =0;
        virtual AbstractRNGPoolElement* getElement(std::size_t index) =0;
    };

    class RandomGeneratorPool: public AbstractRandomGeneratorPool{
    public:
        class RNGPoolElement: public AbstractRNGPoolElement{
        public:
            RNGPoolElement();
            double uniformRealRndValue() override;
            double normalRealRndValue() override;
            MersenneBitSource* getRandomBitSource() override;

        private:
            MersenneBitSource rngGenerator_;
            std::uniform_real_distribution<double> uniformDist_;
            std::normal_distribution<double> normalDist_;
        };

        RandomGeneratorPool();
        RNGPoolElement* getThreadElement();
        RNGPoolElement* getElement(std::size_t index);

    private:
        std::vector<std::unique_ptr<RNGPoolElement>> elements_;
    };


    /**
     * A random distribution which produces random samples with a specific distribution
     */
    class RandomDistribution{
    public:
        virtual double rndValue() =0;
        virtual ~RandomDistribution() = default;
    };

    typedef std::unique_ptr<Core::RandomDistribution> RndDistPtr; ///< pointer type used throughout the project

    /**
     * A uniform distribution, which generates random samples in a specified interval
     */
    class UniformRandomDistribution: public RandomDistribution{
    public:
        UniformRandomDistribution();
        UniformRandomDistribution(double min, double max);
        double rndValue() override;
    private:
        std::uniform_real_distribution<double> internalUniformDist_;
    };

    /**
     * Random distribution which random samples wich are gaussian normal distributed (with mu=0.0 and sigma=1.0)
     */
    class NormalRandomDistribution: public RandomDistribution{
    public:
        NormalRandomDistribution();
        double rndValue() override;
    private:
        std::normal_distribution<double> internalNormalDist_;
    };

    /**
     * A uniform distribution, which generates *non* random test samples in a specified interval.
     * The samples are sequentially chosen from a small set of predetermined test values
     */
    class UniformTestDistribution: public RandomDistribution{
    public:
        UniformTestDistribution() = default;
        UniformTestDistribution(double min, double max);
        double rndValue() override;
    private:
        std::size_t sampleIndex_ = 0;
        double min_ = 0.0;
        double interval_ = 1.0;
    };

    /**
     * Random distribution which *non* random samples wich are gaussian normal distributed (with mu=0.0 and sigma=1.0).
     * The samples are sequentially chosen from a small set of predetermined test values
     */
    class NormalTestDistribution: public RandomDistribution{
    public:
        NormalTestDistribution();
        double rndValue() override;
    private:
        std::size_t sampleIndex_;
    };

    /**
     * Random generator interface: A random generator is able to produce uniform and normal distributed random
     * samples and custom uniform random distributions.
     */
    class AbstractRandomGenerator{
    public:
        virtual ~AbstractRandomGenerator() =default;
        virtual double uniformRealRndValue() =0;
        virtual double normalRealRndValue() =0;
        virtual RndDistPtr getUniformDistribution(double min, double max) =0;
    };

    /**
     * RandomGenerator for production use.
     */
    class PooledRandomGenerator: public AbstractRandomGenerator{
    public:
        PooledRandomGenerator();
        double uniformRealRndValue() override;
        double normalRealRndValue() override;
        std::unique_ptr<RandomDistribution> getUniformDistribution(double min, double max) override;
    };


    /**
     * RandomGenerator for testing: Produces *non* random test values / test distributions
     */
    class TestRandomGenerator: public AbstractRandomGenerator{
    public:
        TestRandomGenerator();
        double uniformRealRndValue() override;
        double normalRealRndValue() override;
        std::unique_ptr<RandomDistribution> getUniformDistribution(double min, double max) override;

    private:
        UniformTestDistribution uniformDistribution_;
        NormalTestDistribution normalDistribution_;
    };

    extern std::unique_ptr<AbstractRandomGenerator> globalRandomGenerator; ///<the global random generator used in the project
    extern std::unique_ptr<AbstractRandomGeneratorPool> globalRandomGeneratorPool; ///< global provider
}

#endif //BTree_randomGenerators
