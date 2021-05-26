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

#include "RS_SimulationConfiguration.hpp"
#include "RS_AbstractReaction.hpp"
#include <sstream>

bool RS::SimulationConfiguration::addSubstance(std::unique_ptr<RS::Substance>& subst) {
    if (substancesNameMap_.count(subst->name()) == 0) {
        substancesNameMap_.insert(
                std::pair<std::string, RS::Substance *>(subst->name(), subst.get())
        );
        if (subst->type() == RS::Substance::substanceType::discrete){
            discreteSubstances_.push_back(subst.get());
        }
        substances_.push_back(std::move(subst));
        return true;
    }
    else{
        return false;
    }
}

RS::Substance* RS::SimulationConfiguration::substance(std::size_t index) {
    return substances_.at(index).get();
}

RS::Substance* RS::SimulationConfiguration::substanceByName(std::string substanceName) {
    try{
        return (substancesNameMap_.at(substanceName));
    }
    catch(std::out_of_range&){
        std::stringstream ss;
        ss << "Substance "<< substanceName<< " is not existing";
        throw (RS::SimulationConfigurationException(ss.str()));
    }
}

std::vector<RS::Substance*> RS::SimulationConfiguration::getAllSubstances() {
    std::vector<RS::Substance*> result = std::vector<RS::Substance*>();
    for(const auto& subst: substances_){
        result.push_back(subst.get());
    }
    return result;
}

std::vector<RS::Substance*> RS::SimulationConfiguration::getAllDiscreteSubstances() {
    return discreteSubstances_;
}


bool RS::SimulationConfiguration::addReaction(std::unique_ptr<RS::AbstractReaction>& reac) {
    reactions_.push_back(std::move(reac));
    return true;
}

RS::AbstractReaction* RS::SimulationConfiguration::reaction(std::size_t index) {
    return reactions_[index].get();
}

std::vector<RS::AbstractReaction*> RS::SimulationConfiguration::getAllReactions() {
    std::vector<RS::AbstractReaction*> result = std::vector<RS::AbstractReaction*>();
    for(const auto& reac: reactions_){
        result.push_back(reac.get());
    }
    return result;
}

void RS::SimulationConfiguration::updateConfiguration() {
    //update all static reaction probabilities in all reactions:
    for(const auto& reac: reactions_){
        reac->updateStaticReactionConcentration();
    }
}