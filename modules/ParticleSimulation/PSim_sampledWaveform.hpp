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
 PSim_swiftWaveform.hpp

 Sampled periodic waveform, for example for SWIFT excitation waveforms
 (with simple input from CSV files)

 ****************************/

#ifndef BTree_swiftWaveform_hpp
#define BTree_swiftWaveform_hpp


#include <vector>
#include <string>

namespace ParticleSimulation{

    /**
     * Sampled waveform (for example for SWIFT simulations)
     * Sampled waveforms are simple one dimensional vectors with voltage values, one sampled value per timestep.
     * Sampled waveforms are provided as simple CSV files, one value per line.
     */
    class SampledWaveform {

    public:
        explicit SampledWaveform(std::string filename);
        [[nodiscard]] bool good() const;
        [[nodiscard]] std::size_t size() const;
        [[nodiscard]] double getValue(std::size_t index) const;
        [[nodiscard]] double getValueLooped(std::size_t index) const;
        [[nodiscard]] double operator[](std::size_t index) const;
        [[nodiscard]] double getInterpolatedValue(double phase) const;

    private:
        std::vector<double> wfTable_;
        std::vector<double> phaseTable_;
        std::size_t size_ = 0;
        bool dataIsGood_ = false;
    };
}

#endif /* BTree_swiftWaveform_hpp */
