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
#include "Core_randomGenerators.hpp"

/**
 * Constructs a RS Simulation with a simulation configuration given by a simulation configuration file
 */
RS::Simulation::Simulation(const std::string& configFileName){
    RS::ConfigFileParser parser = RS::ConfigFileParser();
    std::unique_ptr<RS::SimulationConfiguration> simConf = parser.parseFile(configFileName);
    initFromSimulationConfig_(std::move(simConf));
}

/**
 * Constructs a RS Simulation with a given simulation configuration
 */
RS::Simulation::Simulation(std::unique_ptr<RS::SimulationConfiguration> simConf)
{
    initFromSimulationConfig_(std::move(simConf));
}

RS::SimulationConfiguration* RS::Simulation::simulationConfiguration() const {
    return simConf_.get();
}

bool RS::Simulation::addParticle(RS::ReactiveParticle* particle, index_t index) {
    std::pair<pMap::iterator,bool> insResult = particleMap_.insert(pPair(index,particle));

    if (insResult.second) {
        discreteConcentrations_[particle->getSpecies()]++;
    }
    return insResult.second;
}

void RS::Simulation::removeParticle(index_t index) {
    discreteConcentrations_[particleMap_.at(index)->getSpecies()]--;
    particleMap_.erase(index);
}

RS::ReactiveParticle& RS::Simulation::getParticle(index_t index) const{
    return *particleMap_.at(index);
}

std::map<RS::Substance* const,int> RS::Simulation::discreteConcentrations() const{
    return discreteConcentrations_;
};

/*
 * Gets the number of total reaction events happened in this RS simulation
 */
long RS::Simulation::totalReactionEvents() const{
    return totalReactionEvents_;
}

/*
 * Gets the total number of "ill" events, which are reaction events with a
 * linearized probability => 1
 */
long RS::Simulation::illEvents() const{
    return illEvents_;
}

long RS::Simulation::reactionEvents(RS::AbstractReaction* reaction) const {
    return reactionEvents_.at(reaction);
}

/**
 * Performs a time step with static reaction conditions (same reaction conditions for all particles)
 * @param conditions The reaction conditions for the time step
 * @param dt The time step length
 * @param particleReactedFct An optional function, which is performed for the particles which have reacted
 */
void RS::Simulation::performTimestep(RS::ReactionConditions& conditions, double dt, const particleReactedFctType& particleReactedFct) {

    std::size_t nParticles = particleMap_.size();

    #pragma omp parallel default(none) shared(particleMap_) firstprivate(particleReactedFct, nParticles, conditions, dt)
    {
        reactionMap indMap = indReactDeepCopy_(); // get local copy of independent reaction maps

        if (particleReactedFct != nullptr) {
            #pragma omp for
            for (std::size_t i = 0; i<nParticles; ++i) {
                bool hasReacted = react_(i, conditions, dt, indMap);
                if (hasReacted){
                    RS::ReactiveParticle* particle = particleMap_[i];
                    particleReactedFct(particle);
                }
            }
        }
        else {
            #pragma omp for
            for (std::size_t i = 0; i<nParticles; ++i) {
                react_(i, conditions, dt, indMap);
            }
        }
    }
}

/**
 * Performs a time step with variable reaction conditions for the individual particles, defined by a function.
 * The function takes individual particles and the time at time step begin as parameters.
 *
 * @param conditionFct A function which generates individual reaction conditions for the particles
 * @param dt The time step length
 * @param particleReactedFct An optional function, which is performed for the particles which have reacted
 */
void RS::Simulation::performTimestep(const reactionConditionFctType& conditionFct, double dt, const particleReactedFctType& particleReactedFct) {
    std::size_t nParticles = particleMap_.size();
    double time = this->simulationTime();

    #pragma omp parallel default(none) shared(particleMap_) firstprivate(conditionFct, particleReactedFct, nParticles, time, dt)
    {
        reactionMap indMap = indReactDeepCopy_(); // get local copy of independent reaction maps

        if (particleReactedFct != nullptr) {
            #pragma omp for
            for (std::size_t i = 0; i<nParticles; ++i) {
                RS::ReactiveParticle* particle = particleMap_[i];
                ReactionConditions conditions = conditionFct(particle, time);
                bool hasReacted = react_(i, conditions, dt, indMap);
                if (hasReacted){
                    particleReactedFct(particle);
                }
            }
        }
        else {
            #pragma omp for
            for (std::size_t i = 0; i<nParticles; ++i) {
                RS::ReactiveParticle* particle = particleMap_[i];
                ReactionConditions conditions = conditionFct(particle, time);
                react_(i, conditions, dt, indMap);
            }
        }
    }
}


/**
 * Actually perform a reaction: The particle with index "index" is changed into the product species
 * @param particle the particle to react
 * @param product the new chemical species for this particle
 */
void RS::Simulation::doReaction(RS::AbstractReaction* reaction, RS::ReactiveParticle* particle, RS::Substance* product)
{

    //we have a reaction event: count the reaction event
    #pragma omp atomic
    totalReactionEvents_++;

    #pragma omp atomic
    reactionEvents_[reaction]++;

    // despite the possible illnes of the reaction event: now REACT! => remove actual particle
    // (independent reaction => therefore only on discrete educt and this must be the current particle)
    //this->removeParticle(index);
    //RS::ReactiveParticle* particle = particleMap_.at(index);
    #pragma omp atomic
    discreteConcentrations_[particle->getSpecies()]--;


    //add /append product particles (actually, there is only one possible discrete product, because
    //the number of particles cannot increase (no ion generation in SIMION / RS possible at present)
    //and particle destruction reactions are not yet implemented,
    //an independent reaction has only one discrete educt.
    //(Currently the config file parser enforces only one discrete products)

    //#pragma omp atomic
    particle->setSpecies(product);
    /*RS::ReactiveParticle productParticle = RS::ReactiveParticle(
            product,
            particle.getLocation(),
            particle.getCharge());*/

    //this->addParticle(productParticle,index);
    #pragma omp atomic
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
bool RS::Simulation::react(index_t index, RS::ReactionConditions& conditions, double dt) {
    return react_(index, conditions, dt, reacInd_);
}

bool RS::Simulation::react_(index_t index, RS::ReactionConditions& conditions, double dt, reactionMap &reacInd) {
    RS::ReactiveParticle* particle = particleMap_[index];

    std::vector<AbstractReaction*> &iReactions = reacInd[particle->getSpecies()];
    // shuffle the order of the reaction for every time step to prevent simulation artifacts:
    // std::shuffle(iReactions.begin(), iReactions.end(), *Core::globalRandomGeneratorPool->getThreadRandomSource()->getRandomBitSource());

    for (const auto& reaction: iReactions){ //iterate through all independent reactions

        RS::ReactionEvent reactionEvent = reaction->attemptReaction(conditions, particle, dt);

        if (reactionEvent.reactionHappened){
            if (reactionEvent.reactionProbability > 1.0){ //if the probability is >1 the time step was too long, and the reaction event was ill
                ++illEvents_;
            }
            // despite the possible illness of the reaction event: now REACT!

            RS::Substance* product = reaction->discreteProducts()->begin()->first;

            //#pragma omp critical (do_reaction)
            {
                doReaction(reaction, particle, product);
            }


            //end the independent reactions loop => the actual educt particle does not longer exist
            return true;
        }
    }
    return false;
}

bool RS::Simulation::collisionReact(index_t index, RS::Substance* reactionPartnerSpecies, CollisionConditions& conditions){
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
                RS::Substance* product = reaction->discreteProducts()->begin()->first;
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

int RS::Simulation::timestep() const{
    return nTimesteps_;
}

double RS::Simulation::simulationTime() const{
    return sumTime_;
}

void RS::Simulation::printConcentrations() const{
    std::cout << concentrationString_();
    std::cout << std::endl;
}

void RS::Simulation::printReactionStatistics() const{
    std::cout << reactionStatisticsString_();
    std::cout << std::endl;
}

void RS::Simulation::logConcentrations(std::shared_ptr<spdlog::logger>& logger) const{
    logger->info(concentrationString_());
}

void RS::Simulation::logReactionStatistics(std::shared_ptr<spdlog::logger>& logger) const{
    logger->info(reactionStatisticsString_());
}

std::string RS::Simulation::concentrationString_() const{
    std::stringstream ss;
    ss <<"t= "<<sumTime_<<" ts="<<nTimesteps_<<" | ";
    for(const auto& p: simConf_->getAllDiscreteSubstances()){
        ss << p->name() <<" "<< discreteConcentrations_.at(p)<<" | " ;
    }
    return ss.str();
}

std::string RS::Simulation::reactionStatisticsString_() const{
    std::stringstream ss;
    ss <<"t= "<<sumTime_<<" ts="<<nTimesteps_<<" | \n";
    for(const auto& reaction: simConf_->getAllReactions()){
        ss << reaction->getLabel() <<": "<< reactionEvents_.at(reaction)<<" |"<<std::endl;
    }
    return ss.str();
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

void RS::Simulation::initFromSimulationConfig_(std::unique_ptr<RS::SimulationConfiguration> simConf) {
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
            RS::Substance* particleSubstance = nullptr;
            RS::Substance* backgroundSubstance = nullptr;

            for(const auto& substPair: *reac->educts()){
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
            auto discreteEduct = reac->discreteEducts()->begin();
            reacInd_.at(discreteEduct->first).push_back(reac);
            //staticProbabilities_.at(discreteEduct->first).push_back(reac->staticReactionConcentration());
        }
        else{
            //the reaction is dependent on more than one discrete educt: add it to the reaction maps of all educts
            for(const auto& discreteEduct: *reac->discreteEducts()){
                reacDep_.at(discreteEduct.first).push_back(reac);
            }
        }
    }
}

RS::Simulation::reactionMap RS::Simulation::indReactDeepCopy_() {
    reactionMap result;

    auto it = reacInd_.begin();
    while(it != reacInd_.end())
    {
        result[it->first] = it->second;
        ++it;
    }

    return result;
}
