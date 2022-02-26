/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2022 - Physical and Theoretical Chemistry /
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

 FMM3D_fmmSolver_C_interface.h

 Coulombic Particle / Particle (space charge) solver based on FMM3D library:
 C++ interface/C function wrapper to separate C++ compiler from incompatible C code in FMM3D

 ****************************/

#ifndef IDSIMF_FMM3D_FMMSOLVER_C_INTERFACE_H
#define IDSIMF_FMM3D_FMMSOLVER_C_INTERFACE_H

void lfmm3d_s_c_g_wrapper(double *eps, int *nsource,
                   double *source, double *charge, double *pot, double *grad, int *ier);

#endif //IDSIMF_FMM3D_FMMSOLVER_C_INTERFACE_H
