/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2022 - Physical and Theoretical Chemistry /
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
 PSim_sampledFunction.hpp

 Simple one dimensional sampled function: Defines a one dimensional function
 given by a set of discrete samples, with simple CSV input

 ****************************/
#ifndef IDSIMF_PSIM_SAMPLEDFUNCTION_HPP
#define IDSIMF_PSIM_SAMPLEDFUNCTION_HPP

#include <vector>
#include <string>

namespace ParticleSimulation{

    /**
     * Sampled function (for example for ion mobility scalings)
     * Sampled functions are simple one dimensional vectors with discrete function samples, given as pairs of
     * independent and dependent values
     * Sampled functions are provided as simple CSV files, with independent / dependent values per line,
     * separated by a delimiter
     */
    class SampledFunction {

    public:
        explicit SampledFunction(std::string filename);
        [[nodiscard]] bool good() const;
        [[nodiscard]] std::size_t size() const;
        [[nodiscard]] double getFunctionValue(std::size_t index) const;
        [[nodiscard]] double getIndependentValue(std::size_t index) const;
        [[nodiscard]] std::pair<double,double> operator[](std::size_t index) const;
        [[nodiscard]] double getInterpolatedValue(double independentValue) const;

    private:
        const std::string delimiter_ = ";"; ///<A delimiter for the columns in the function sample input files
        std::vector<double> functionValsTable_;
        std::vector<double> independentValsTable_;
        std::size_t size_ = 0;
        bool dataIsGood_ = false;
    };
}

#endif //IDSIMF_PSIM_SAMPLEDFUNCTION_HPP
