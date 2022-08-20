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
 ****************************/
#include "PSim_sampledFunction.hpp"
#include <fstream>
#include <cmath>
#include <algorithm>


ParticleSimulation::SampledFunction::SampledFunction(std::string filename) {
    //check if file with filename is accessible, then read data:
    dataIsGood_ = false;

    std::ifstream fIn(filename.c_str());
    if (fIn.good()){

        // fill tables with values from csv:
        size_t pos = 0;
        std::string line;
        std::string leftColumn;

        while(std::getline(fIn, line)) {
            if (line[0] != '#') {
                //find a token terminated by delimiter_, extract it and delete it from the read line
                if ((pos = line.find(delimiter_)) != std::string::npos) {
                    leftColumn = line.substr(0, pos);
                    line.erase(0, pos + delimiter_.length());
                    independentValsTable_.push_back(std::strtod(leftColumn.c_str(), nullptr));
                    functionValsTable_.push_back(std::strtod(line.c_str(), nullptr));
                }
                else {
                    throw std::invalid_argument("Line without delimiter between columns found in function sample file");
                }
            }
        }

        if (independentValsTable_.size() != functionValsTable_.size()){
            throw std::invalid_argument("Different lengths of function samples and independent sample vectors found");
        }
        size_ = independentValsTable_.size();
        dataIsGood_ = true;
    }
}

/**
 * Checks if the input file was read correctly and function data is ready
 */
bool ParticleSimulation::SampledFunction::good() const{
    return dataIsGood_;
}

/**
 * Returns the size of the function values vector
 */
std::size_t ParticleSimulation::SampledFunction::size() const{
    return size_;
}

/**
 * Gets a value from the function values vector
 * @param index a numeric index in the function values vector of the function
 */
double ParticleSimulation::SampledFunction::getFunctionValue(std::size_t index) const{
    return functionValsTable_.at(index);
}

/**
 * Gets a value from the independent values ("x" axis) vector
 * @param index a numeric index in the independent values vector of the function
 */
double ParticleSimulation::SampledFunction::getIndependentValue(std::size_t index) const {
    return independentValsTable_.at(index);
}

/**
 * Sampled functions can be used like arrays
 */
std::pair<double,double> ParticleSimulation::SampledFunction::operator[](std::size_t index) const{
    return {getIndependentValue(index), getFunctionValue(index)};
}

/**
 * Gets linearly interpolated function values
 *
 * @param independentValue The independent value ("x"-Value) to interpolate the function value (dependent value) for
 */
double ParticleSimulation::SampledFunction::getInterpolatedValue(double independentValue) const{

    auto upper = std::upper_bound(independentValsTable_.begin(), independentValsTable_.end(), independentValue);
    if (upper == independentValsTable_.begin() || upper == independentValsTable_.end() ){
        throw std::invalid_argument("Value to interpolate for not within independent value range of the sampled function");
    }

    // calculate nearest two function samples
    size_t upperIndex = static_cast<std::size_t>((--upper) - independentValsTable_.begin());

    double indepLower = independentValsTable_.at(upperIndex);
    double indepUpper = independentValsTable_.at(upperIndex+1);

    double funcLower = functionValsTable_.at(upperIndex);
    double funcUpper = functionValsTable_.at(upperIndex+1);

    double weight = (independentValue-indepLower) /(indepUpper - indepLower);
    double funcInterpolated = funcLower*(1-weight) + funcUpper*weight;

    return funcInterpolated;
}