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
#include "BTree_abstractNode.hpp"

#include "BTree_abstractNode.hpp"
#include "BTree_particle.hpp"
#include <sstream>
#include <iostream>

/**
 * Constructs a tree node
 * @param min the lower corner of the node (spatial xlo,ylo,zlo corner)
 * @param max the upper corner of the node (spatial xhi,yhi,zhi corner)
 * @param parent the parent node, the new node will be a subnode of that node (can be nullptr for root nodes)
 */
BTree::AbstractNode::AbstractNode(Core::Vector min, Core::Vector max)
        :
        min_(min),
        max_(max),
        center_(min+(max-min)/2)
{}

/**
 * Copy constructor
 * @param that a node to copy
 */
BTree::AbstractNode::AbstractNode(const BTree::AbstractNode& that):
        charge_(that.charge_),
        centerOfCharge_(that.centerOfCharge_),
        numP_(that.numP_),
        min_(that.min_),
        max_(that.max_),
        center_(that.center_),
        particle_(that.particle_)
{}

// Static methods:

/**
 * Return total number of existing tree nodes
 */
int BTree::AbstractNode::getNumberOfNodes(){
    return nNodes_;
}

/**
 * Computes the electric field on a test charge located at r1 from a charge located at r2
 * @param r1 the position of the small test charge (in m)
 * @param r2 the position of the charge "charge2" (in m)
 * @param charge2 the charge of "charge2" (in coulomb)
 * @return The electric field at the position r1
 */
Core::Vector BTree::AbstractNode::calculateElectricField(const Core::Vector &r1, const Core::Vector &r2, double charge2){
    double d = (r1-r2).magnitude();
    if (d>0){
        return ( (r1-r2)*(charge2/(Core::ELECTRIC_CONSTANT*(d*d*d))) ); //d*d*d is faster than pow(d,3.0)
    }
    else{
        return (Core::Vector(0.0,0.0,0.0));
    }
}



// Accessors:

/**
 * Gets the number of particles in the node including the suboctant nodes and thus all particles in the
 * spatial extend of the node
 */
std::size_t BTree::AbstractNode::getNumberOfParticles() const{
    return numP_;
}

/**
 * Gets a pointer to the particle in the node
 */
BTree::Particle* BTree::AbstractNode::getParticle() const{
    return particle_;
}

/**
 * Returns the spatial center of the node
 */
Core::Vector BTree::AbstractNode::getCenter() const{
    return center_;
}

/**
 * Returns the "lower" corner of the node (corner with xLo, yLo, zLo coordinates)
 */
Core::Vector BTree::AbstractNode::getMin() const{
    return min_;
}

/**
 * Returns the "upper" corner of the  node (corner with xHi, yHi, zHi coordinates)
 */
Core::Vector BTree::AbstractNode::getMax() const{
    return max_;
}

/**
 * Returns the total charge of the node (including the contributions of the subnodes)
 */
double BTree::AbstractNode::getCharge() const{
    return charge_;
}

/**
 * Returns the geometric center of charge of the node (including the contribution of the subnodes)
 */
Core::Vector BTree::AbstractNode::getCenterOfCharge() const{
    return centerOfCharge_;
}

/**
 * Returns the sub-octant of the current node where a given location falls into
 * @param location a location to get the  sub-octant for
 * @return the sub octant in which given location falls into
 */
BTree::AbstractNode::Octant BTree::AbstractNode::getOctant(Core::Vector location) const{
    if (location.z()<=center_.z()){
        if(location.x()<=center_.x()){
            if(location.y()<=center_.y()){
                return SWB;
            }
            else{
                return NWB;
            }
        }
        else{
            if(location.y()>=center_.y()){
                return NEB;
            }
            else {
                return SEB;
            }
        }
    }
    else{
        if(location.x()<=center_.x()){
            if(location.y()<=center_.y()){
                return SWT;
            }
            else{
                return NWT;
            }
        }
        else{
            if(location.y()>=center_.y()){
                return NET;
            }
            else{
                return SET;
            }
        }
    }
}