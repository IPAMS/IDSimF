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
 ****************************/

#include "FMM3D_fmmSolver.hpp"
extern "C" {
    #include "FMM3D_fmmSolver_C_interface.h"
}
#include "Core_randomGenerators.hpp"
#include <vector>
#include <iostream>

Core::Vector FMM3D::FMMSolver::computeEFieldFromSpaceCharge(Core::Particle& particle) {

    computeChargeDistribution();
    return{1.0, 0.0, 0.0};
}


void FMM3D::FMMSolver::computeChargeDistribution() {
    Core::RandomSource* rndSource = Core::globalRandomGeneratorPool->getThreadRandomSource();

    int ns=200000;
    int ier=0;

    std::vector<double> sources(3*ns);
    std::vector<double> charges(ns);
    std::vector<double> pot(ns);
    std::vector<double> grad(3*ns);

    // initialize arrays
    for(int i=0; i<ns; ++i)
    {
        sources[3*i] = pow(rndSource->uniformRealRndValue() ,2);
        sources[3*i+1] = pow(rndSource->uniformRealRndValue() ,2);
        sources[3*i+2] = pow(rndSource->uniformRealRndValue() ,2);

        charges[i] = 1.0;

    }

    double eps = 0.5e-6;

// call the fmm routine
    lfmm3d_s_c_g_wrapper(&eps, &ns, sources.data(), charges.data(), pot.data(), grad.data(), &ier);

    std::cout << "FMM Result" << pot[0]<<" "<<grad[0]<<" "<<grad[1]<<" "<<grad[2];
}