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
 RS_ConcentrationFileWriter.hpp

 This implements a file writer for tabular concentration result data from an RS simulation

 ****************************/

#ifndef RS_ConcentrationFileWriter_hpp
#define RS_ConcentrationFileWriter_hpp

#include <fstream>
#include "RS_SimulationConfiguration.hpp"
#include "RS_Simulation.hpp"
#include <vector>

namespace RS {
    class ConcentrationFileWriter {

    public:
        explicit ConcentrationFileWriter (std::string transientFilename);

        void initFile(SimulationConfiguration* simConf);
        void closeFile();
        void writeTimestep(Simulation& sim);
        void writeReactionStatistics(Simulation& sim);

    private:
        std::ofstream transientFile_;
        static const std::string partSeparatorString;
    };
}

#endif //RS_ConcentrationFileWriter_hpp
