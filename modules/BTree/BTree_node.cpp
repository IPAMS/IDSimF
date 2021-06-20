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

#include "BTree_node.hpp"
#include "BTree_particle.hpp"
#include <iostream>


/**
 * Constructs a tree node
 * @param min the lower corner of the node (spatial xlo,ylo,zlo corner)
 * @param max the upper corner of the node (spatial xhi,yhi,zhi corner)
 * @param parent the parent node, the new node will be a subnode of that node (can be nullptr for root nodes)
 */
BTree::Node::Node(Core::Vector min, Core::Vector max, BTree::Node* parent):
        GenericBaseNode<Node>(min,max,parent)
{}

// Public member methods:
/**
 * Computes the electric field on a given particle from the charged particles in a tree
 * with this node as root node
 * @param targetP a particle to calculate the total electric field for
 * @return the electric force on the particle resulting from the particles in the tree with this node as root
 */
Core::Vector BTree::Node::computeElectricFieldFromTree(BTree::Particle &targetP){
    if (numP_ == 1){
        Core::Vector efield= calculateElectricField(
                targetP.getLocation(),
                particle_->getLocation(),
                particle_->getCharge());
        return(efield);
    }
    else{
        double r = (targetP.getLocation()-centerOfCharge_).magnitude();
        double d = max_.x() - min_.x();
        
        if (r > 0 && d/r < BTree::Node::theta){
            Core::Vector efield= calculateElectricField(
                    targetP.getLocation(),
                    centerOfCharge_,
                    charge_);
            return(efield);
            
        }
        else{
            Core::Vector efield= Core::Vector(0.0,0.0,0.0);
            for (auto & octNode : octNodes_){
                if(octNode != nullptr){
                    efield = efield +octNode->computeElectricFieldFromTree(targetP);
                }
            }
            return(efield);
        }
    }
}


/**
 * Tests if all tree nodes are within the spatial limits of their parent nodes
 *
 * @throws logic_error if an inconsistent node is found
 */
void BTree::Node::testSpatialTreeIntegrity(){
    for (auto & octNode : octNodes_){
        if(octNode != nullptr){
            octNode->testSpatialTreeIntegrity();
        }
    }

    if (this->parent_ != nullptr){
        
        if ( (this->min_.x() < this->parent_->min_.x()) ||
             (this->min_.y() < this->parent_->min_.y()) ||
             (this->min_.z() < this->parent_->min_.z()) ||
             (this->max_.x() > this->parent_->max_.x()) ||
             (this->max_.y() > this->parent_->max_.y()) ||
             (this->max_.z() > this->parent_->max_.z())
            ){
            std::stringstream ss;
            ss << " >>>>> inconsistent node found this-min:"<<this->min_<<" max "<<this->max_ << " parent-min:"<<this->parent_->min_<<" max "<<this->parent_->max_<<std::endl;
            throw (std::logic_error(ss.str()));
        }
    }
}


/**
 * Tests if the particles in the nodes of the subtree with this node as root are consistent, regarding the
 * spatial positions.
 *
 * @throws logic_error if a particle is found which is out of the spatioal bounds of its containing node
 */
void BTree::Node::testNodeParticleIntegrity(){
    
    if (this->particle_ != nullptr){
        Core::Vector pLoc = this->particle_->getLocation();
        if (
            pLoc.x() < this->min_.x() ||
            pLoc.y() < this->min_.y() ||
            pLoc.z() < this->min_.z() ||
            pLoc.x() > this->max_.x() ||
            pLoc.y() > this->max_.y() ||
            pLoc.z() > this->max_.z()
            ){
            
            std::stringstream ss;
            ss << "Node with illegal particle found : p: "<<this->particle_<<std::endl;
            ss << this->toString();
            throw (std::logic_error(ss.str()));
        }
    }
    
    for (auto & octNode : octNodes_){
        if(octNode != nullptr){
            octNode->testNodeParticleIntegrity();
        }
    }
}

/**
 * Checks if a tree node is in a subtree with this node as root
 *
 * @param nodeToFind pointer to the node to find in the subtree
 * @param debug prints information about the found node to cout if true
 * @return true if the given node is found in the subtree
 */
bool BTree::Node::isNodeInSubtree(const BTree::Node *nodeToFind, bool debug) const{
    bool result = false;
    if (this == nodeToFind){
        if (debug) {
            std::cout << " >>>>> node found this-min:" << this->min_ << " max " << this->max_ << " nPart" << this->numP_
                      << " " << std::endl;
        }
        return true;
    }

    for (auto & octNode : octNodes_){
        if(octNode != nullptr){
            if (octNode->isNodeInSubtree(nodeToFind, debug)){
                result = true;
            };
        }
    }
    return result;
}

/**
 * Checks if a particle is in a subtree with this node as root
 *
 * @param particle pointer to a particle to find in the subtree
 * @param debug prints informatio about the found particle if true
 * @return true if the particle is found in the subtree
 */
bool BTree::Node::isParticleInSubtree(const BTree::Particle *particle, bool debug) const{

    bool result = false;
    if (this->particle_ == particle){
        if(debug) {
            std::cout << " >>>> particle found <<<<" << std::endl;
            std::cout << this->toString();
        }
        return true;
    }
    
    
    for (auto & octNode : octNodes_){
        if(octNode != nullptr){
            if (octNode->isParticleInSubtree(particle, debug) == true){
                result = true;
            };
        }
    }
    return result;
}