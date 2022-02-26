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

#include "BTree_parallelTree.hpp"
#include "Core_utils.hpp"
#include <iostream>

/**
 Constructor: Constructs a new tree
 
 \param min the lower spatial boundary of the tree domain
 \param min the upper spatial boundary of the tree domain
 */
BTree::ParallelTree::ParallelTree(Core::Vector min, Core::Vector max):
        root_(std::make_unique<BTree::ParallelNode>(min,max,nullptr)) {
    root_->initAsRoot();
    iVec_ = std::make_unique<treeParticlePtrList>();
    iMap_ = std::make_unique<std::unordered_map<std::size_t, treeParticlePtrList::const_iterator>>();
}

/* copy operator missing
 */
//BTRee::Tree

/** 
  Get the tree root node
 */
BTree::ParallelNode* BTree::ParallelTree::getRoot() const{
    return(root_.get());
}

/**
 Get a linear particle linked list of the particles in the tree
 */
BTree::treeParticlePtrList* BTree::ParallelTree::getParticleList() const{
    return(iVec_.get());
}

/**
 Gets the number of particles in the tree
 */
std::size_t BTree::ParallelTree::getNumberOfParticles() const{
    return(root_->getNumberOfParticles());
}

/**
 * Inits the internal data structures after a structural change on the
 * tree structure and returns the total number of tree nodes.
 *
 * @return the total number of nodes in this tree
 */
std::size_t BTree::ParallelTree::init() {
    return updateNodes(1);
}

/**
 * Counts and returns the number of nodes on the individual tree levels and
 */
std::vector<std::size_t> BTree::ParallelTree::countNodesOnLevels()
{
    std::vector<std::size_t> numOfNodesOnLevels(nTreeLevels_);
    root_->countNodesOnLevel(1, numOfNodesOnLevels);
    return numOfNodesOnLevels;
}

/**
 * Computes the electric field acting on a given particle from the charged particles in this tree
 * @param particle A particle to calculate the total electric field for
 * @return The electric field on the particle resulting from the particles in the tree
 */
Core::Vector BTree::ParallelTree::getEFieldFromSpaceCharge(Core::Particle &particle){

    Core::Vector efield= Core::Vector(0.0,0.0,0.0);
    Core::Vector loc=particle.getLocation();

    std::vector<BTree::ParallelNode*>::size_type cur=0; // index for the current position in the working array

    // construct a linear node list to store the nodes to process in:
    std::vector<BTree::ParallelNode*> nodesToProcess;
    nodesToProcess.reserve(nodesReserveSize_);
    // add the tree root to the process list to init the calculation process:
    nodesToProcess.push_back(root_.get());

    while(cur < nodesToProcess.size())
    {
        //process the current node in the node process list
        BTree::ParallelNode* currentNode = nodesToProcess.at(cur);

        if(currentNode->numP_ == 1){ // if node has only one particle: calculate force directly
            efield=efield+root_->calculateElectricField(loc, currentNode->particle_->wrappedParticle->getLocation(),
                    currentNode->particle_->wrappedParticle->getCharge());
        }
        else { // if more particles: process the node
            // get squared distance to the charge center:
            double r2= (loc-nodesToProcess[cur]->centerOfCharge_).magnitudeSquared();

            // if distance is large enough: Use node as approximate pseudo particle (use charge center and charge weight of the node)
            if(currentNode->edgeLengthSquaredNormalized_ < r2){
                efield=efield+root_->calculateElectricField(loc, currentNode->centerOfCharge_, currentNode->charge_);
            }
            else { // if distiance is too small: Decent a level in the tree and add the child nodes to the process list
                for(auto & octNode : currentNode->octNodes_) {
                    if(octNode != nullptr){
                        nodesToProcess.push_back(octNode);
                    }
                }
            }
        }
        // Increase the current node index to process the next node. If no nodes were added to the node process list,
        // the current index is at the end of the process list and the loop is terminated.
        cur++;
    }
    return(efield);
}

/**
 Insert particle into the tree
 
 \param particle the particle to insert
 \param ext_index an external index number for the particle / numerical particle id (most likely from simion)
 */
void BTree::ParallelTree::insertParticle(Core::Particle &particle, std::size_t ext_index){

    auto treeParticle = std::make_unique<BTree::TreeParticle>(&particle);
    root_->insertParticle(treeParticle.get());
    iVec_->push_front(std::move(treeParticle));
    iMap_->insert({ext_index, iVec_->cbegin()});
}

/**
 Removes a particle with a given particle index / particle id and its hosting leaf node
 from the tree
 
 \param ext_index the external numerical particle id
 */
void BTree::ParallelTree::removeParticle(std::size_t ext_index){
    
    auto iter =(*iMap_)[ext_index];
    BTree::TreeParticle* particle = iter->get();//[int_index];
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
BTree::TreeParticle* BTree::ParallelTree::getParticle(std::size_t ext_index) const{
    auto iter =(*iMap_)[ext_index];
    BTree::TreeParticle* particle = iter->get();
    return (particle);
}

/**
 * Upates the location of a particle in this tree.
 * @param extIndex Index / ID of the Particle to modify
 * @param newLocation Location to set for the selected particle
 * @param numNodesChanged an integer reference to count the number of structurally changed nodes
 */
void BTree::ParallelTree::updateParticleLocation(std::size_t extIndex, Core::Vector newLocation, int* numNodesChanged){

    BTree::TreeParticle* particle = getParticle(extIndex);
    BTree::AbstractNode* pNode = particle->getHostNode();

    if ( (   newLocation.x() <= pNode->getMin().x() ||
             newLocation.y() <= pNode->getMin().y() ||
             newLocation.z() <= pNode->getMin().z() )
         ||
         (   newLocation.x() >= pNode->getMax().x() ||
             newLocation.y() >= pNode->getMax().y() ||
             newLocation.z() >= pNode->getMax().z() ) )
    {
        (*numNodesChanged)++;
        Core::Particle* wrappedParticle = particle->wrappedParticle;
        removeParticle(extIndex);
        wrappedParticle->setLocation(newLocation);
        insertParticle(*wrappedParticle, extIndex);
    }
    else{
        particle->wrappedParticle->setLocation(newLocation);
    }
}

/**
 * Updates the internal state of the nodes in this tree. The update is performed
 * serialized and parallelized.
 *
 * @param ver number of nodes which have changed in terms of tree structure
 * (0 means that the structure of the tree has not changed and the serialized
 * data structures has not to be updated)
 * @return the total number of nodes in the tree
 */
std::size_t BTree::ParallelTree::updateNodes(int ver) {
    if(ver!=0)
    {
        //if a structural change in the tree structure has happened: We have to update serialized
        //data structures
        nTreeLevels_ = getTreeDepth_();
        nodeStartIndicesOnLevels_.resize(nTreeLevels_, 0);
        nodesOnLevels_ = countNodesOnLevels();
        numberOfNodesTotal_=0;
        for(std::size_t i=0; i<nTreeLevels_; ++i)
        {
            numberOfNodesTotal_+=nodesOnLevels_[i];
        }
        nodesReserveSize_ = numberOfNodesTotal_ / 4;

        updateLevelStartIndices_();
        serializeNodes_();
    }

    updateNodeChargeState_();
    return numberOfNodesTotal_;
}


void BTree::ParallelTree::printParticles() const{
    int i =0;
    for (auto&& particle : *iVec_) {
        
        std::cout<<"particle "<<i <<" " << particle->wrappedParticle->getLocation() <<" "<<particle->getHostNode()->getMin() << " "<< particle->getHostNode()->getMax()<<std::endl;
        i++;
    }
    
    std::cout <<" keys:";
    for(std::pair<const int, BTree::treeParticlePtrList::const_iterator> kv : *iMap_) {
        std::cout<< kv.first <<" | ";
    }
    std::cout <<std::endl;
}


/**
 * Get number of tree levels in this tree
 */
std::size_t BTree::ParallelTree::getTreeDepth_() const {
    return root_->maximumRecursionDepth();
}

/**
 * Updates the the internal vector of indices to the first nodes on the individual levels of the tree
 */
void BTree::ParallelTree::updateLevelStartIndices_()
{
    nodeStartIndicesOnLevels_[0]=0;
    for(std::size_t i=1; i<nTreeLevels_; ++i)
    {
        nodeStartIndicesOnLevels_.at(i)=nodeStartIndicesOnLevels_.at(i-1)+nodesOnLevels_.at(i-1);
    }
}

/**
 * Retrieves and serializes the links to tree nodes into nodesSerialized_ vector
 */
void BTree::ParallelTree::serializeNodes_()
{
    nodesSerialized_.resize(numberOfNodesTotal_, nullptr);

    //copy the start indices vector to a working copy since
    //the serialization method modifies the passed indices vector
    std::vector<std::size_t> currentStartIndices = nodeStartIndicesOnLevels_;
    root_->serializeIntoVector(nodesSerialized_, 0, currentStartIndices);
}

void BTree::ParallelTree::updateNodeChargeState_()
{
    for(std::size_t i=nTreeLevels_; i>=1; i--)
    {
        // get boundaries of the current tree level in the linearized vector to process.
        // the deepest tree level has the end of the linearized vector as upper limit,
        // the other have the begin of the next level as upper limit
        std::size_t levelUpperBound = i==nTreeLevels_ ? numberOfNodesTotal_ : nodeStartIndicesOnLevels_.at(i);
        std::size_t levelLowerBound = nodeStartIndicesOnLevels_.at(i-1);

        #pragma omp parallel for
        for(std::size_t j = levelUpperBound; j > levelLowerBound; j--)
        {
            ParallelNode* currentNode = nodesSerialized_.at(j-1);

            if(currentNode->numP_==1) {
                // If the current node has only one particle: Update node parameters with parameters from particle
                currentNode->centerOfCharge_ = currentNode->particle_->wrappedParticle->getLocation();
                currentNode->charge_ = currentNode->particle_->wrappedParticle->getCharge();
            }
            else {
                // If the current node represents multiple particles, and has therefore sub nodes,
                // update node parameters from child nodes.
                // Since we start updating with the lowest level, we know that the sub nodes are already
                // updated.
                currentNode->charge_ = 0.0;
                currentNode->centerOfCharge_ = Core::Vector(0.0,0.0,0.0);
                for (auto & octNode : currentNode->octNodes_){
                    if(octNode != nullptr){
                        currentNode->charge_ += octNode->getCharge();
                        currentNode->centerOfCharge_ = currentNode->centerOfCharge_+
                                                        (octNode->centerOfCharge_ *
                                                         octNode->charge_);
                    }
                }

                if (Core::isDoubleUnequal(currentNode->charge_, 0.0)){
                    currentNode->centerOfCharge_ = currentNode->centerOfCharge_ / currentNode->charge_;
                }
                else {
                    currentNode->centerOfCharge_ = currentNode->center_;
                }
            }
        }
    }
}