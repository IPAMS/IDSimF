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

    /**
 * A random distribution which produces random samples with a specific distribution
 */
    class RandomDistribution{
    public:
        virtual double rndValue() =0;
        virtual ~RandomDistribution() = default;
    };

    using RndDistPtr= std::unique_ptr<Core::RandomDistribution>; ///< pointer type used throughout the project

    /**
     * A uniform distribution, which generates random samples in a specified interval
     */
    class UniformRandomDistribution: public RandomDistribution{
    public:
        UniformRandomDistribution(double min, double max, RandomBitSource<rndBit_type>* randomSource);
        double rndValue() override;
    private:
        RandomBitSource<rndBit_type>* randomSource_;
        std::uniform_real_distribution<double> internalUniformDist_;
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

    class AbstractRandomGeneratorPool{
    public:
        virtual ~AbstractRandomGeneratorPool() =default;
        virtual RndDistPtr getUniformDistribution(double min, double max) =0;
        virtual AbstractRNGPoolElement* getThreadElement() =0;
        virtual AbstractRNGPoolElement* getElement(std::size_t index) =0;
    };

    class RandomGeneratorPool: public AbstractRandomGeneratorPool{
    public:
        class RNGPoolElement: public AbstractRNGPoolElement{
        public:
            RNGPoolElement() = default;
            double uniformRealRndValue() override;
            double normalRealRndValue() override;
            MersenneBitSource* getRandomBitSource() override;

        private:
            MersenneBitSource rngGenerator_;
            std::uniform_real_distribution<double> uniformDist_;
            std::normal_distribution<double> normalDist_;
        };

        RandomGeneratorPool();
        RndDistPtr getUniformDistribution(double min, double max) override;
        RNGPoolElement* getThreadElement() override;
        RNGPoolElement* getElement(std::size_t index) override;

    private:
        std::vector<std::unique_ptr<RNGPoolElement>> elements_;
    };

    class TestRandomGeneratorPool: public AbstractRandomGeneratorPool{
    public:
        class TestPoolElement: public AbstractRNGPoolElement{
        public:
            TestPoolElement() = default;
            double uniformRealRndValue() override;
            double normalRealRndValue() override;
            TestBitSource* getRandomBitSource() override;

        private:
            TestBitSource rngGenerator_;
            UniformTestDistribution uniformDist_;
            NormalTestDistribution normalDist_;
        };

        TestRandomGeneratorPool() = default;
        RndDistPtr getUniformDistribution(double min, double max) override;
        TestPoolElement* getThreadElement() override;
        TestPoolElement* getElement(std::size_t index) override;

    private:
        TestPoolElement element_;

    };

    extern std::unique_ptr<AbstractRandomGeneratorPool> globalRandomGeneratorPool; ///< global provider
}

#endif //BTree_randomGenerators
