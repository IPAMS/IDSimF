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
 ****************************/

#include "appUtils_stopwatch.hpp"

void AppUtils::Stopwatch::start() {
    beginCpu_ = std::clock();
    beginWall_ = std::chrono::system_clock::now();
}

void AppUtils::Stopwatch::stop() {
    clock_t endCpu = std::clock();
    elapsedSecsCPU_ += double(endCpu - beginCpu_) / CLOCKS_PER_SEC;
    elapsedWall_ += (std::chrono::system_clock::now() - beginWall_);
}

double AppUtils::Stopwatch::elapsedSecondsCPU() {
    return elapsedSecsCPU_;
}

double AppUtils::Stopwatch::elapsedSecondsWall() {
    return elapsedWall_.count();
}

