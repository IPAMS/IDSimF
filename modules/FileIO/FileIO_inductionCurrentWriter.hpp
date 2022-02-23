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
 FileIO_inductionCurrentWriter.hpp

 Simulation result file writer which writes induced current on arbitrary detection electrodes, defined
 by weight fields given as potential arrays

 ****************************/

#ifndef Particle_simulation_qit_fft_writer
#define Particle_simulation_qit_fft_writer

#include <fstream>
#include <vector>
#include <memory>

//forward declare own classes:
namespace Core{
    class Particle;
}

namespace BTree{
    class Vector;
}

namespace ParticleSimulation{
    class SimionPotentialArray;
}
    /**
     *  Simulation result file writer which writes the displacement current induced by a simulated ion cloud
     *  on arbitrary detection electrodes.
     *  The induction effect on the individual electrodes is defined by individual potential arrays, which contain
     *  the field which would result with the individual detection electrode on a normalized voltage (1V) and
     *  all other electrodes in the geometry on ground.
     *  The effect of the individual electrodes can be scaled by a scaling factor, to allow bipolar or more complex
     *  detection circuits of the detection electrodes.
     */

namespace FileIO{
    class InductionCurrentWriter{

    public:
        InductionCurrentWriter(std::vector<Core::Particle *> particles, std::string transientFilename,
                               const std::vector<ParticleSimulation::SimionPotentialArray*> &weightFields,
                               std::vector<double> weightFactors,
                               double scale_mm_per_gu);
        ~InductionCurrentWriter();

        void writeTimestep(double time);

    private:
        double scale_mm_per_gu_ = 0.0;
        std::unique_ptr<std::ofstream> transientFile_;
        std::vector<Core::Particle*> particles_;
        std::vector<std::pair<ParticleSimulation::SimionPotentialArray*, double>> weightFields_;
    };
}


#endif //Particle_simulation_qit_fft_writer
