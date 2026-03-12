/***************************
 Ion Dynamics Simulation Framework (IDSimF)

 Copyright 2026 - Physical and Theoretical Chemistry /
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
#include "CollisionModel_HDF5MDTrajectoryWriter.hpp"

CollisionModel::HDF5MDTrajectoryWriter::HDF5MDTrajectoryWriter(std::string hdf5Filename) {
    h5f_ = std::make_unique<FileIO::HDF5File>(hdf5Filename, FileIO::HDF5File::WRITE_ONLY);
}

void CollisionModel::HDF5MDTrajectoryWriter::initNewTrajectory(
        CollisionModel::Molecule &molecule, CollisionModel::Molecule &bgMolecule) {
    nTrajectories_++;
    std::size_t nAtomsMolecule = molecule.getAtoms().size();
    std::size_t nAtomsBG = bgMolecule.getAtoms().size();

    currentDS_ = h5f_->initTableDataset("MD_trajectories", "trajectory"+std::to_string(nTrajectories_),
        (nAtomsMolecule+nAtomsBG)*3 + 2);

    std::vector<std::string> columnNames={"time", "dt"};
    std::vector<double> masses;

    for (std::size_t i=0; i<nAtomsMolecule; i++) {
        columnNames.push_back("mol_a"+std::to_string(i)+"_pos_x");
        columnNames.push_back("mol_a"+std::to_string(i)+"_pos_y");
        columnNames.push_back("mol_a"+std::to_string(i)+"_pos_z");
        masses.push_back(molecule.getAtoms()[i]->getMass() / Core::AMU_TO_KG);
    }

    for (std::size_t i=0; i<nAtomsBG; i++) {
        columnNames.push_back("bg_a"+std::to_string(i)+"_pos_x");
        columnNames.push_back("bg_a"+std::to_string(i)+"_pos_y");
        columnNames.push_back("bg_a"+std::to_string(i)+"_pos_z");
        masses.push_back(bgMolecule.getAtoms()[i]->getMass() / Core::AMU_TO_KG);
    }
    h5f_->writeDatasetAttribute<std::string>(currentDS_,"column_names", columnNames);
    h5f_->writeDatasetAttribute<double>(currentDS_,"masses", masses);

    std::vector<std::size_t> nAtoms = {nAtomsMolecule, nAtomsBG};
    h5f_->writeDatasetAttribute<std::size_t>(currentDS_,"number_of_atoms", nAtoms);
}

void CollisionModel::HDF5MDTrajectoryWriter::writeTrajectorySample(
    double time, double dt, Molecule& molecule, Molecule& bgMolecule) {

    std::vector<double> row = {time, dt };
    auto atomsMolecule = molecule.getAtoms();
    for (const auto& atom : atomsMolecule) {
        // Core::Vector atomPos = atom->getRelativePosition() + molecule.getComPos();
        Core::Vector atomPos = molecule.genWorldFramePosition(atom);
        row.push_back(atomPos.x());
        row.push_back(atomPos.y());
        row.push_back(atomPos.z());
    }

    auto atomsBg = bgMolecule.getAtoms();
    for (const auto& atom : atomsBg) {
        Core::Vector atomPos = bgMolecule.genWorldFramePosition(atom);
        row.push_back(atomPos.x());
        row.push_back(atomPos.y());
        row.push_back(atomPos.z());
    }

    h5f_->writeRowToTableDataset(currentDS_, row);
}
