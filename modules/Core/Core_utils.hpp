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
 Core_utils.hpp

 Utility functions

 ****************************/
#ifndef IDSIMF_CORE_UTILS_HPP
#define IDSIMF_CORE_UTILS_HPP

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
namespace Core {
    inline bool isDoubleUnequal(double lhs, double rhs){
        return lhs != rhs;
    }

    inline bool isDoubleEqual(double lhs, double rhs){
        return lhs == rhs;
    }
}
#pragma GCC diagnostic pop

#endif //IDSIMF_CORE_UTILS_HPP
