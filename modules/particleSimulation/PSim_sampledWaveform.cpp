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

#include "PSim_sampledWaveform.hpp"
#include <fstream>
#include <cmath>

/**
 * Constructor: Creates a sampled waveform from a given file
 * @param filename the name of the file to read
 */
ParticleSimulation::SampledWaveform::SampledWaveform(std::string filename){
    //check if file with filename is accessible, then read data:
    dataIsGood_ = false;

    std::ifstream fIn(filename.c_str());
    if (fIn.good()){

        // fill wave table with values from csv:
        std::string line;
        while(std::getline(fIn, line)) {
            if (line[0] != '#') {
                wfTable_.push_back(std::strtod(line.c_str(), nullptr));
            }
        }
        size_ = wfTable_.size();

        // fill phases table:
        double phaseIncrement = 1.0 / size_;
        for(std::size_t i=0; i<size_; ++i){
            phaseTable_.push_back(i*phaseIncrement);
        }

        dataIsGood_ = true;
    }
}

/**
 * Checks if the input file was read correctly and waveform data is ready
 * @return true if the data was read correctly
 */
bool ParticleSimulation::SampledWaveform::good() const{
    return dataIsGood_;
}

/**
 * Returns the size of the waveform vector
 *
 * @return the number of elements in the waveform vector
 */
std::size_t ParticleSimulation::SampledWaveform::size() const{
    return size_;
}

/**
 * Gets a value from the waveform data
 * @param index an numeric index in the data vector of the waveform
 * @return Voltage value of the timestep specified by the index
 */
double ParticleSimulation::SampledWaveform::getValue(std::size_t index) const{
    return wfTable_.at(index);
}

/**
 * Gets a value from the waveform data in a looped way: If the index is larger
 * than the size of the waveform vector, the index is wrapped around
 * Example: Index 12 in a 10 element waveform gets index 1 (not 2 because of zero indexing)
 * @param index an numeric index in the data vector of the waveform
 * @return Voltage value of the timestep specified by the index
 */
double ParticleSimulation::SampledWaveform::getValueLooped(std::size_t index) const{
    return this->getValue(index % this->size_);
}

/**
 * Sampled waveforms can be used like arrays
 * @param index an numeric index in the data vector of the sampled waveform
 * @return Voltage value specified by index
 */
double ParticleSimulation::SampledWaveform::operator[](std::size_t index) const{
    return this->getValue(index);
}

/**
 * Gets linearly interpolated waveform values
 *
 * @param phase The waveform phase to get an interpolated value for (phase is considered not to be in radians, thus
 * the phase values have to be between [0 .. 1])
 * @return The linearly interpolated value for the given phase
 */
double ParticleSimulation::SampledWaveform::getInterpolatedValue(double phase) {
    // calculate nearest two values

    std::size_t iLower = static_cast<std::size_t>(std::floor(phase * size_));
    std::size_t iHigher = iLower+1;


    //"normal" interpolation between samples:
    if (iHigher < size_) {
        double phase0 = phaseTable_[iLower];
        double phase1 = phaseTable_[iHigher];
        double val0 = wfTable_[iLower];
        double val1 = wfTable_[iHigher];

        // calculate interpolated value
        return val0 + (phase - phase0) * ((val1-val0) / (phase1-phase0));
    }
    //edge case: Interpolation between last and first sample of the waveform
    else if (iHigher == size_){
        double phase0 = phaseTable_[iLower];
        double phase1 = 1.0;
        double val0 = wfTable_[iLower];
        double val1 = wfTable_[0];

        // calculate interpolated value
        return val0 + (phase - phase0) * ((val1-val0) / (phase1-phase0));
    }
    else if (iHigher > size_){
        throw (std::invalid_argument("Illegal phase (out of valid data range)"));
    }
}