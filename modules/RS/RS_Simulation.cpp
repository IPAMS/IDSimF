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

#include "RS_Simulation.hpp"

RS::Simulation::Simulation(std::string configFileName){
    RS::ConfigFileParser parser = RS::ConfigFileParser();
    RS::Simulation(parser.getTestConfigWaterClusters());
}

RS::Simulation::Simulation(std::unique_ptr<RS::SimulationConfiguration> simConf):
        nTimesteps_(0),
        sumTime_(0.0),
        totalReactionEvents_(0),
        illEvents_(0)
{

    simConf_ = std::move(simConf);
    substances_ = simConf_->getAllSubstances();
    reactions_ = simConf_->getAllReactions();

    for(const auto& subst: substances_){
        //std::cout << subst <<std::endl;
        reacInd_[subst] = std::vector<AbstractReaction*>();
        reacDep_[subst] = std::vector<AbstractReaction*>();
        if (subst->type() == RS::Substance::substanceType::discrete){
            discreteConcentrations_[subst] = 0;
        }
    }

    for(const auto& reac: reactions_){

        //prepare reaction event counter:
        reactionEvents_.insert({reac, 0});

        if (reac->isCollisionReaction()){
            //the reaction is a collision based reaction: there is only one discrete educt and one non discrete educt
            RS::Substance* particleSubstance;
            RS::Substance* backgroundSubstance;

            for(const auto& substPair: reac->educts()){
                RS::Substance* subst = substPair.first;
                int factor = substPair.second;
                if (factor == 1){
                    if (subst->type() == RS::Substance::discrete){
                        particleSubstance = subst;
                    }else{
                        backgroundSubstance = subst;
                    }
                }
            }
            //add reaction to collision reaction table
            std::pair<RS::Substance*,RS::Substance*> eductPair(particleSubstance,backgroundSubstance);
            if (reacCollision_.count(eductPair) == 0){
                reacCollision_[eductPair] = std::vector<AbstractReaction*>();
            }
            reacCollision_.at(eductPair).push_back(reac);
            //throw (RS::ConfigurationFileException("Illegal reaction in configuration file"));
        }
        else if (reac->isIndependent()){
            //the reaction is independent: there is only one discrete educt
            auto discreteEduct = reac->discreteEducts().begin();
            reacInd_.at(discreteEduct->first).push_back(reac);
            //staticProbabilities_.at(discreteEduct->first).push_back(reac->staticReactionConcentration());
        }
        else{
            //the reaction is dependent on more than one discrete educt: add it to the reaction maps of all educts
            for(const auto& discreteEduct: reac->discreteEducts()){
                reacDep_.at(discreteEduct.first).push_back(reac);
            }
        }
    }
}

RS::SimulationConfiguration* RS::Simulation::simulationConfiguration() {
    return simConf_.get();
}

bool RS::Simulation::addParticle(RS::ReactiveParticle* particle,int index) {
    std::pair<pMap::iterator,bool> insResult = particleMap_.insert(pPair(index,particle));

    if (insResult.second) {
        discreteConcentrations_[particle->getSpecies()]++;
    }
    return insResult.second;
}

void RS::Simulation::removeParticle(int index) {
    discreteConcentrations_[particleMap_.at(index)->getSpecies()]--;
    particleMap_.erase(index);
}

RS::ReactiveParticle& RS::Simulation::getParticle(int index) {
    return *particleMap_.at(index);
}

std::map<RS::Substance* const,int> RS::Simulation::discreteConcentrations(){
    return discreteConcentrations_;
};

/*
 * Gets the number of total reaction events happened in this RS simulation
 */
long RS::Simulation::totalReactionEvents(){
    return totalReactionEvents_;
}

/*
 * Gets the total number of "ill" events, which are reaction events with a
 * linearized probability => 1
 */
long RS::Simulation::illEvents(){
    return illEvents_;
}

long RS::Simulation::reactionEvents(RS::AbstractReaction* reaction) const {
    return reactionEvents_.at(reaction);
}

/**
 * Actually perform a reaction: The particle with index "index" is changed into the product species
 * @param particle the particle to react
 * @param product the new chemical species for this particle
 */
void RS::Simulation::doReaction(RS::AbstractReaction* reaction, RS::ReactiveParticle* particle, RS::Substance* product)
{

    //we have a reaction event: count the reaction event
    totalReactionEvents_++;
    reactionEvents_[reaction]++;

    // despite the possible illnes of the reaction event: now REACT! => remove actual particle
    // (independent reaction => therefore only on discrete educt and this must be the current particle)
    //this->removeParticle(index);
    //RS::ReactiveParticle* particle = particleMap_.at(index);
    discreteConcentrations_[particle->getSpecies()]--;


    //add /append product particles (actually, there is only one possible discrete product, because
    //the number of particles cannot increase (no ion generation in SIMION / RS possible at present)
    //and particle destruction reactions are not yet implemented,
    //an independent reaction has only one discrete educt.
    //(Currently the config file parser enforces only one discrete products)

    particle->setSpecies(product);
    /*RS::ReactiveParticle productParticle = RS::ReactiveParticle(
            product,
            particle.getLocation(),
            particle.getCharge());*/

    //this->addParticle(productParticle,index);
    discreteConcentrations_[product]++;
}

/**
 * Let a particle react: The independent reactions of that particle are tested if they occur,
 * if one reaction occurs, that reaction is performed
 * @param index the index of the particle to react
 * @param conditions the parameters present while the reaction occurs
 * @param dt the time step length
 * @return true if a reaction had occurred
 */
bool RS::Simulation::react(int index, RS::ReactionConditions& conditions, double dt) {
    RS::ReactiveParticle* particle = particleMap_.at(index);
    std::vector<AbstractReaction*> iReactions = reacInd_.at(particle->getSpecies());

    for(const auto& reaction: iReactions){ //iterate through all independent reactions

        RS::ReactionEvent reactionEvent = reaction->attemptReaction(conditions, particle, dt);

        if (reactionEvent.reactionHappened){
            if (reactionEvent.reactionProbability > 1.0){ //if the probability is >1 the time step was too long, and the reaction event was ill
                ++illEvents_;
            }
            // despite the possible illnes of the reaction event: now REACT!

            RS::Substance* product = reaction->discreteProducts().begin()->first;
            doReaction(reaction, particle, product);

            //end the independent reactions loop => the actual educt particle does not longer exist
            return true;
        }
    }
    return false;
}

bool RS::Simulation::collisionReact(int index, RS::Substance* reactionPartnerSpecies, CollisionConditions& conditions){
    //step 1: get reactions for this particle and this reaction partner
    RS::ReactiveParticle* particle = particleMap_.at(index);
    RS::Substance* particleSpecies = particle->getSpecies();
    std::pair<RS::Substance*,RS::Substance*> eductPair(particleSpecies,reactionPartnerSpecies);

    //step 2: Iterate through all possible reactions and check if they happen...
    if (reacCollision_.count(eductPair) > 0){
        std::vector<AbstractReaction*> reactions= reacCollision_.at(eductPair);

        //iterate through all collision based reactions for this pair of particle and collision partner:
        for (const auto &reaction: reactions) {
            RS::ReactionEvent reactionEvent = reaction->attemptReaction(conditions, particle);

            //if probability is >= 1.0 => do reaction
            if (reactionEvent.reactionHappened){
                RS::Substance* product = reaction->discreteProducts().begin()->first;
                doReaction(reaction, particle, product);
                return true;
            }
        }
    }
    return false;
}

void RS::Simulation::advanceTimestep(double dt) {
    sumTime_ += dt;
    nTimesteps_ ++;
}

int RS::Simulation::timestep() {
    return nTimesteps_;
}

double RS::Simulation::simulationTime() {
    return sumTime_;
}

void RS::Simulation::printConcentrations() {
    std::cout <<"t= "<<sumTime_<<" ts="<<nTimesteps_<<" | ";
    for(const auto& p: simConf_->getAllDiscreteSubstances()){
        std::cout << p->name() <<" "<< discreteConcentrations_.at(p)<<" | " ;
    }
    std::cout << std::endl;
}

void RS::Simulation::printReactionStatistics() {
    std::cout <<"t= "<<sumTime_<<" ts="<<nTimesteps_<<" | \n";
    for(const auto& reaction: simConf_->getAllReactions()){
        std::cout << reaction->getLabel() <<": "<< reactionEvents_.at(reaction)<<" |"<<std::endl;
    }
    std::cout << std::endl;
}


std::ostream& operator<<(std::ostream& os, const RS::Simulation& sim)
{
    os <<"RS_Simulation, substances: "<<std::endl;
    for(const auto& subst: sim.substances_){
        os << *subst <<std::endl;
    }

    os << "reactions:"<<std::endl;
    for(const auto& reac: sim.reactions_){
        os << *reac <<std::endl;
    }

    os << "independent reactions:"<<std::endl;
    for(const auto& sub_reac: sim.reacInd_){
        os << sub_reac.first <<std::endl;
        for(const auto& reac: sub_reac.second){
            os << *reac << std::endl;
        }
    }

    os << "dependent reactions:"<<std::endl;
    for(const auto& sub_reac: sim.reacDep_){
        os << sub_reac.first <<std::endl;
        for(const auto& reac: sub_reac.second){
            os << *reac << std::endl;
        }
    }

    os <<"=============== particles in simulation =================="<<std::endl;
    for(const auto& particlePair: sim.particleMap_){
        os << particlePair.first << " | " << particlePair.second <<std::endl;
    }

    os <<"=============== concentrations in simulation =================="<<std::endl;
    for(const auto& p: sim.discreteConcentrations_){
        os << *p.first << " | " << p.second <<std::endl;
    }

    return os;
}