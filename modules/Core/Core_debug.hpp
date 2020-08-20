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
 Core_debug.hpp

 Basic global debugging switch / debugging facilities

 ****************************/
#ifndef IDSIMF_CORE_DEBUG_HPP
#define IDSIMF_CORE_DEBUG_HPP

namespace Core {
    inline bool debug = false;
    inline bool safetyGuards = true;  //< global switch for additional runtime safety guards
}

#endif //IDSIMF_CORE_DEBUG_HPP
