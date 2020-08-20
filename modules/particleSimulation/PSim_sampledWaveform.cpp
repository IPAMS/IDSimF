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

/**
 * Constructor: Creates a sampled waveform from a given file
 * @param filename the name of the file to read
 */
ParticleSimulation::SampledWaveform::SampledWaveform(std::string filename){
    //check if file with filename is accessible, then read data:
    dataIsGood_ = false;

    std::ifstream fIn(filename.c_str());
    if (fIn.good()){
        std::string line;
        while(std::getline(fIn, line)) {
            if (line[0] != '#') {
                wfTable_.push_back(std::strtod(line.c_str(), nullptr));
            }
        }
        size_ = wfTable_.size();
        dataIsGood_ = true;
    }
}

/**
 * Checks if the input file was read correctly and waveform data is ready
 * @return true if the data was read correctly
 */
bool ParticleSimulation::SampledWaveform::good(){
    return dataIsGood_;
}

/**
 * Returns the size of the waveform vector
 *
 * @return the number of elements in the waveform vector
 */
std::size_t ParticleSimulation::SampledWaveform::size(){
    return size_;
}

/**
 * Gets a value from the waveform data
 * @param index an numeric index in the data vector of the waveform
 * @return Voltage value of the timestep specified by the index
 */
double ParticleSimulation::SampledWaveform::getValue(std::size_t index){
    return wfTable_.at(index);
}

/**
 * Gets a value from the waveform data in a looped way: If the index is larger
 * than the size of the waveform vector, the index is wrapped around
 * Example: Index 12 in a 10 element waveform gets index 1 (not 2 because of zero indexing)
 * @param index an numeric index in the data vector of the waveform
 * @return Voltage value of the timestep specified by the index
 */
double ParticleSimulation::SampledWaveform::getValueLooped(std::size_t index){
    return this->getValue(index % this->size_);
}

/**
 * Sampled waveforms can be used like arrays
 * @param index an numeric index in the data vector of the sampled waveform
 * @return Voltage value specified by index
 */
double ParticleSimulation::SampledWaveform::operator[](std::size_t index){
    return this->getValue(index);
}
