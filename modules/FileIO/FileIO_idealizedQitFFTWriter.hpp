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
 FileIO_idealizedQitFFTWriter.hpp

 FFT Spectrum writer for idealized QIT Simulations

 Simulation result file writer which writes the approximate electric field on
 the cap electrodes in an idealized QIT simulation

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

namespace FileIO{
    /**
     *  Simulation result file writer which writes the displacement current induced by a simulated ion cloud
     *  on the cap electrodes in an QIT simulation to a file.
     *
     *  The displacement current is approximated as the average displacement of the individual charged particles.
     */
    class IdealizedQitFFTWriter{

    public:
        IdealizedQitFFTWriter(std::vector<Core::Particle*> particles, std::string transientFilename);
        ~IdealizedQitFFTWriter();

        void writeTimestep(double time);
        void writeTimestepMassResolved(double time);

    private:
        std::unique_ptr<std::ofstream> transientFile_;
        std::vector<Core::Particle*> particles_;
    };
}


#endif //Particle_simulation_qit_fft_writer
