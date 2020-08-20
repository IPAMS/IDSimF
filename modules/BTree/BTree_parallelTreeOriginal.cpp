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

#include "BTree_parallelTreeOriginal.hpp"
#include <iostream>

/**
 Constructor: Constructs a new tree
 
 \param min the lower spatial boundary of the tree domain
 \param min the upper spatial boundary of the tree domain
 */
BTree::ParallelTreeOriginal::ParallelTreeOriginal(Core::Vector min, Core::Vector max){
    root_ = std::make_unique<BTree::ParallelNodeOriginal>(min,max,nullptr);
    root_->initAsRoot();
    
    iVec_ = std::make_unique<std::list<Particle*>>();
    iMap_ = std::make_unique<std::unordered_map<int,std::list<Particle*>::const_iterator>>();
}

/* copy operator missing
 */
//BTRee::Tree

/** 
  Get the tree root
 
 \returns the root Node of the tree
 */
BTree::ParallelNodeOriginal* BTree::ParallelTreeOriginal::getRoot(){
    return(root_.get());
}

/**
 Get a linear particle list of the particles in the tree
 
 \returns a linearized linked list of the particles in the tree
 */
std::list<BTree::Particle*>* BTree::ParallelTreeOriginal::getParticleList(){
    return(iVec_.get());
}

/**
 Get the number of particles in the tree
 
 \returns the number of particles in the tree
 */
int BTree::ParallelTreeOriginal::getNumberOfParticles(){
    return(root_->getNumberOfParticles());
}


/**
 Insert particle into the tree
 
 \param particle the particle to insert
 \param ext_index an external index number for the particle / numerical particle id (most likely from simion)
 */
void BTree::ParallelTreeOriginal::insertParticle(BTree::Particle &particle, int ext_index){
    
    root_->insertParticle(&particle);
    iVec_->push_front(&particle);
    iMap_->insert({ext_index, iVec_->cbegin()});
}

/**
 Removes a particle with a given particle index / particle id and its hosting leaf node
 from the tree
 
 \param ext_index the external numerical particle id
 */
void BTree::ParallelTreeOriginal::removeParticle(int ext_index){
    
    std::list<Particle*>::const_iterator iter =(*iMap_)[ext_index];
    BTree::Particle* particle = *iter;//[int_index];
    BTree::AbstractNode* pHostNode = particle->getHostNode();
    pHostNode->removeMyselfFromTree();
    if (!pHostNode->isRoot()){
        delete(pHostNode);
    }
    
    iMap_->erase(ext_index);
    iVec_->erase(iter);
}

/** 
 Retrieves a particle from the tree by its external index
 
 \param ext_index the external particle index
 \returns the retrieved particle 
 */
BTree::Particle* BTree::ParallelTreeOriginal::getParticle(int ext_index){
    std::list<Particle*>::const_iterator iter =(*iMap_)[ext_index];
    BTree::Particle* particle = *iter;
    return (particle);
}

//FIXME: remove diagnose
void BTree::ParallelTreeOriginal::updateParticleLocationnew(int ext_index, Core::Vector newLocation,int* ver){

    BTree::Particle* particle = this->getParticle(ext_index);
    BTree::AbstractNode* pNode = particle->getHostNode();

    //FIXME: reformat / comment
    if ( (
                 newLocation.x() <= pNode->getMin().x() ||
                 newLocation.y() <= pNode->getMin().y() ||
                 newLocation.z() <= pNode->getMin().z()
                )
                 ||
                (
                 newLocation.x() >= pNode->getMax().x() ||
                 newLocation.y() >= pNode->getMax().y() ||
                 newLocation.z() >= pNode->getMax().z()
                ) 
    )
    {
        (*ver)++;
        this->removeParticle(ext_index);
        particle->setLocation(newLocation);
        this->insertParticle(*particle, ext_index);
    }
    else{
        particle->setLocation(newLocation);
    }


}


void BTree::ParallelTreeOriginal::bestimme_max_rek_tiefe(int* n)
{
    *n=root_->max_rek_tiefe(root_.get());
}

void BTree::ParallelTreeOriginal::bestimme_anzahl_Knoten_auf_ebene(int ebene, std::vector<int> &zaehler, int n)
{
    root_->anz_Knoten_auf_der_Ebene(root_.get(),ebene,zaehler,n);
}

void BTree::ParallelTreeOriginal::bestimme_Feld_begin(std::vector<int> &begin, std::vector<int> &zaehler, int n)
{
    root_->bestimme_Anfang(begin,zaehler,n);
}

void BTree::ParallelTreeOriginal::einfuegen(std::vector<BTree::ParallelNodeOriginal*> &Nodes, int ebene, std::vector<int> &zaehler, std::vector<int> &begin, int n)
{
    root_->insertNodes(root_.get(),Nodes,ebene,zaehler,begin,n);
}

void BTree::ParallelTreeOriginal::Knotenup(std::vector<BTree::ParallelNodeOriginal*> &Nodes, std::vector<int> &begin, int n, int summe)
{
    root_->update(Nodes,begin,n,summe);
}

void BTree::ParallelTreeOriginal::updateS(int ext_index, Core::Vector newLocation)
{
    BTree::Particle* particle = this->getParticle(ext_index);
    BTree::AbstractNode* pNode = particle->getHostNode();
    
    if (particle->getLocation() == newLocation){
        return;
    }
    else{
        if ( !(
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
        //std::cout<< "i="<<ext_index<<" setLocS "<< newLocation<<std::endl;
            particle->setLocation(newLocation);
        }
        else {
            std::cout << "THIS SHOULD NOT HAPPEN"<<std::endl;
        }
    }
}


//todo: Update particle location method (with particle location check etc.)
//todo: Test particle location update
//commit and implement simple SIMION interface and perform simple simulations

/**
 Computes / initalizes the charge distribution in the tree: The charge in the 
 individual tree nodes is calculated from the particles in the tree
 */
void BTree::ParallelTreeOriginal::computeChargeDistribution(){
    root_->computeChargeDistributionRecursive();
}

void BTree::ParallelTreeOriginal::printParticles(){
    int i =0;
    for (BTree::Particle* particle : *iVec_) {
        
        std::cout<<"particle "<<i <<" " << particle->getLocation() <<" "<<particle->getHostNode()->getMin() << " "<< particle->getHostNode()->getMax()<<std::endl;
        i++;
    }
    
    std::cout <<" keys:";
    for(std::pair<const int, std::list<Particle*>::const_iterator> kv : *iMap_) {
        std::cout<< kv.first <<" | ";
    }
    std::cout <<std::endl;
}

/**
 Computes the electric field from the particles in the tree on a given 
 test particle
 
 \param particle a testparticle on which the electric force from the tree acts
 \returns the electric force on that particle
 */
Core::Vector BTree::ParallelTreeOriginal::computeEFieldFromTree(BTree::Particle &particle,std::vector<BTree::ParallelNodeOriginal*> &MyNodes){
    return root_->computeEFieldFromTree(particle,MyNodes);
}
