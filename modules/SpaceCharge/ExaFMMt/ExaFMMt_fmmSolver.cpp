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

#include "ExaFMMt_fmmSolver.hpp"
#include "build_tree.h"
#include "build_list.h"
#include "laplace.h"

void ExaFMMt::FMMSolver::computeChargeDistribution(){

    // step 1: Prepare sources and targets
    std::size_t nParticles = getNumberOfParticles();

    exafmm_t::Bodies<exafmm_t::real_t> sources(nParticles);
    exafmm_t::Bodies<exafmm_t::real_t> targets(nParticles);

    std::size_t i = 0;
    Core::Vector particlePos;
    std::vector<SpaceCharge::particleListEntry*> pListEntriesPtrs;
    for (SpaceCharge::particleListEntry& pListEntry : *iVec_) {
        pListEntriesPtrs.emplace_back(&pListEntry);
        sources[i].ibody = i;
        targets[i].ibody = i;
        particlePos = pListEntry.particle -> getLocation();
        sources[i].X[0] = particlePos.x();
        sources[i].X[1] = particlePos.y();
        sources[i].X[2] = particlePos.z();

        targets[i].X[0] = particlePos.x();
        targets[i].X[1] = particlePos.y();
        targets[i].X[2] = particlePos.z();

        sources[i].q = pListEntry.particle -> getCharge();
        ++i;
    }

    // step 2: Create an Fmm instance for Laplace kernel.
    int P = 8;         // expansion order
    int ncrit = 400;   // max number of bodies per leaf
    exafmm_t::LaplaceFmm fmm(P, ncrit);

    // step 3: Build and balance the octree.
    exafmm_t::NodePtrs<exafmm_t::real_t> leafs, nonleafs;
    exafmm_t::Nodes<exafmm_t::real_t> nodes;

    exafmm_t::get_bounds(sources, targets, fmm.x0, fmm.r0);
    nodes = exafmm_t::build_tree(sources, targets, leafs, nonleafs, fmm);


    // step 4: Build lists and pre-compute invariant matrices.
    exafmm_t::init_rel_coord();
    exafmm_t::build_list(nodes, fmm);
    fmm.M2L_setup(nonleafs);
    fmm.precompute();

    // step 5: Use FMM to evaluate potential
    fmm.upward_pass(nodes, leafs, false);
    fmm.downward_pass(nodes, leafs, false);

    //#pragma omp parallel for default(none) shared(leafs, pListEntriesPtrs)
    for (size_t k=0; k<leafs.size(); ++k) {
        exafmm_t::Node<exafmm_t::real_t>* leaf = leafs[k];
        std::vector<int> & itrgs = leaf->itrgs;
        for (size_t j=0; j<itrgs.size(); ++j) {
            std::cout << itrgs[j] <<" ";
            pListEntriesPtrs[itrgs[j]]->potential = leaf->trg_value[4*j+0];
            pListEntriesPtrs[itrgs[j]]->gradient = {
                    leaf->trg_value[4*j+1],
                    leaf->trg_value[4*j+2],
                    leaf->trg_value[4*j+3]
                    };
        }
    }
}