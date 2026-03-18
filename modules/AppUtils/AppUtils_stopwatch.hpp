/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2021 - Physical and Theoretical Chemistry /
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
 appUtils_stopwatch.hpp

 Simulation timing / runtime measurement

 ****************************/

#ifndef IDSIMF_APPUTILS_STOPWATCH_HPP
#define IDSIMF_APPUTILS_STOPWATCH_HPP

#include <ctime>
#include <chrono>

namespace AppUtils{

    using sysClTimePt = std::chrono::time_point<std::chrono::system_clock>;

    class Stopwatch {
    public:
        void start();
        void stop();
        double elapsedSecondsCPU();
        double elapsedSecondsWall();

    private:
        clock_t beginCpu_;
        double elapsedSecsCPU_ = 0;
        std::chrono::duration<double> elapsedWall_ = std::chrono::duration<double>::zero();

        sysClTimePt beginWall_;

    };
}

#endif //IDSIMF_APPUTILS_STOPWATCH_HPP
