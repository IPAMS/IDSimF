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
 ****************************/

#include "RS_ConcentrationFileWriter.hpp"
const std::string RS::ConcentrationFileWriter::partSeparatorString= "----";

RS::ConcentrationFileWriter::ConcentrationFileWriter(std::string transientFilename) {
    transientFile_ = std::ofstream();
    transientFile_.open(transientFilename);
}

void RS::ConcentrationFileWriter::initFile(RS::SimulationConfiguration *simConf) {
    //write header:
    transientFile_<<"RS C++ result"<<std::endl;
    transientFile_<<"Timestep ; Time ; ";
    std::vector<RS::Substance*> discreteSubst = simConf->getAllDiscreteSubstances();
    for (auto subst: discreteSubst){
        transientFile_<<" "<<subst->name()<<" ;";
    }
    transientFile_<<std::endl;
}

void RS::ConcentrationFileWriter::closeFile() {
    transientFile_.flush();
    transientFile_.close();
}

void RS::ConcentrationFileWriter::writeTimestep(Simulation& sim) {
    std::vector<Substance*> discreteSubst = sim.simulationConfiguration()->getAllDiscreteSubstances();
    std::map<Substance* const,int> concs = sim.discreteConcentrations();

    transientFile_<<sim.timestep()<<" ; "<<sim.simulationTime() <<" ; ";
    for (auto subst: discreteSubst){
        transientFile_<<" "<<concs.at(subst)<<" ;";
    }
    transientFile_ << std::endl;
}

void RS::ConcentrationFileWriter::writeReactionStatistics(Simulation& sim) {
    std::vector<AbstractReaction*> reactions = sim.simulationConfiguration()->getAllReactions();

    transientFile_<<partSeparatorString << std::endl << "Reaction Events | ";
    for (auto reaction: reactions){
        transientFile_<< reaction->getLabel()<<": "<<sim.reactionEvents(reaction)<<" | " ;
    }
    transientFile_ << std::endl;
}