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

#include "BTree_tree.hpp"
#include "Core_particle.hpp"
#include <iostream>

/**
 * Constructs a new tree
 * @param min the lower spatial boundary of the tree domain
 * @param max the upper spatial boundary of the tree domain
 */
BTree::Tree::Tree(Core::Vector min, Core::Vector max):
        root_(std::make_unique<BTree::Node>(min,max,nullptr)){
    root_->initAsRoot();
    iVec_ = std::make_unique<treeParticlePtrList>();
    iMap_ = std::make_unique<std::unordered_map<std::size_t, treeParticlePtrList::const_iterator>>();
}

/**
 * Gets the tree root
 */
BTree::Node* BTree::Tree::getRoot() const {
    return(root_.get());
}

/**
 * Gets a linear particle list of the particles in the tree
 * @return a linearized linked list of the particles in the tree
 */
BTree::treeParticlePtrList* BTree::Tree::getParticleList(){
    return(iVec_.get());
}

/**
 * Get the number of particles in the tree
 *
 * @returns the number of particles in the tree
 */
std::size_t BTree::Tree::getNumberOfParticles() const{
    return(root_->getNumberOfParticles());
}

/**
 * Insert an unwrapped particle into the tree
 *
 * @param particle the particle to insert
 * @param ext_index an external index number for the particle / numerical particle id (most likely from a SIMION simulation)
 */
void BTree::Tree::insertParticle(Core::Particle &particle, size_t ext_index){
    insertParticle_(particle, ext_index);
}

/**
 * Insert an unwrapped particle into the tree
 *
 * @param particle the particle to insert
 * @param ext_index an external index number for the particle / numerical particle id (most likely from a SIMION simulation)
 */
BTree::TreeParticle* BTree::Tree::insertParticle_(Core::Particle &particle, size_t ext_index){

    //create new wrapped particle for tree:
    auto treeParticle = std::make_unique<BTree::TreeParticle>(&particle);
    BTree::TreeParticle* treeParticlePtr = treeParticle.get();
    root_->insertParticle(treeParticle.get());
    iVec_->push_front(std::move(treeParticle));
    iMap_->insert({ext_index, iVec_->cbegin()});
    return treeParticlePtr;
}

/**
 * Removes a particle with a given particle index / particle id and its hosting leaf node
 * from the tree
 *
 * @param ext_index the external numerical particle index
 */
void BTree::Tree::removeParticle(size_t ext_index){
    
    auto iter =(*iMap_)[ext_index];
    BTree::TreeParticle* particle = iter->get();
    BTree::AbstractNode* pHostNode = particle->getHostNode();
    pHostNode->removeMyselfFromTree();
    if (!pHostNode->isRoot()){
        delete(pHostNode);
    }
    
    iMap_->erase(ext_index);
    iVec_->erase(iter);
}

/** 
 * Retrieves a particle from the tree by its external index
 *
 * @param ext_index the external particle index
 * @returns the retrieved particle
 */
BTree::TreeParticle* BTree::Tree::getParticle(size_t ext_index) const{
    auto iter =(*iMap_)[ext_index];
    BTree::TreeParticle* particle = iter->get();
    return (particle);
}

/**
 * Updates the location of a particle in the tree
 * (the position of particles managed by the tree should not be updated directly via the particle,
 * this will invalidate the tree)
 *
 * @param ext_index the external index of the particle to update
 * @param newLocation a new location of that particle
 */
void BTree::Tree::updateParticleLocation(size_t ext_index, Core::Vector newLocation){
    //test if location has changed
    //if yes: test if particle is still in the node

    BTree::TreeParticle* particle = this->getParticle(ext_index);
    BTree::AbstractNode* pNode = particle->getHostNode();
    
    if (particle->get()->getLocation() == newLocation){
        return;
    }
    else if ( !(
                 newLocation.x() <= pNode->getMin().x() ||
                 newLocation.y() <= pNode->getMin().y() ||
                 newLocation.z() <= pNode->getMin().z()
                )
                 &&
                !(
                 newLocation.x() >= pNode->getMax().x() ||
                 newLocation.y() >= pNode->getMax().y() ||
                 newLocation.z() >= pNode->getMax().z()
                )  )
    {
        particle->get()->setLocation(newLocation);
        pNode->updateSelf();
        pNode->updateParents();
    }
    else {
        // we have to reinsert
        Core::Particle* wrappedParticle = particle->get();
        this->removeParticle(ext_index); //this destroys the reference to particle (the TreeParticle)
        wrappedParticle->setLocation(newLocation);
        BTree::TreeParticle* newTreeParticle = this->insertParticle_(*wrappedParticle, ext_index);
        newTreeParticle->getHostNode()->updateSelf();
        newTreeParticle->getHostNode()->updateParents();
    }
}

/**
 * Computes / initalizes the charge distribution in the tree: The charge in the
 * individual tree nodes is calculated from the particles in the tree
 */
void BTree::Tree::computeChargeDistribution(){
    root_->computeChargeDistributionRecursive();
}

/**
 * Computes the electric field from the particles in the tree on a given
 * test particle
 *
 * @param particle a testparticle on which the electric force from the tree acts
 * @returns the electric force on that particle
 */
Core::Vector BTree::Tree::computeEFieldFromTree(Core::Particle &particle){
    return root_->computeElectricFieldFromTree(particle);
}

/**
 * Prints the particles in the tree to cout
 */
void BTree::Tree::printParticles() const{
    int i =0;
    for (auto&& particle : *iVec_ ) {

        std::cout<<"particle "<<i <<" " << particle->get()->getLocation() <<" "<<particle->getHostNode()->getMin() << " "<< particle->getHostNode()->getMax()<<std::endl;
        i++;
    }

    std::cout <<" keys:";
    for(std::pair<const int, treeParticlePtrList::const_iterator> kv : *iMap_) {
        std::cout<< kv.first <<" | ";
    }
    std::cout <<std::endl;
}