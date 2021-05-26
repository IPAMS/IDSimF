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

#include "BTree_parallelNode.hpp"
#include <math.h>
#include <sstream>
#include <iostream>

/**
 * Constructs a tree node
 * @param min the lower corner of the node (spatial xlo,ylo,zlo corner)
 * @param max the upper corner of the node (spatial xhi,yhi,zhi corner)
 * @param parent the parent node, the new node will be a subnode of that node (can be nullptr for root nodes)
 */
BTree::ParallelNode::ParallelNode(Core::Vector min, Core::Vector max, BTree::ParallelNode* parent):
        GenericBaseNode<ParallelNode>(min,max,parent)
{}


/**
 * Copy constructor
 */
BTree::ParallelNode::ParallelNode(const BTree::ParallelNode& that):
        GenericBaseNode<ParallelNode>(that)
{}


/**
 * Asssignment operator
 */
BTree::ParallelNode& BTree::ParallelNode::operator=(const BTree::ParallelNode& that){

    GenericBaseNode::operator=(that);
    return *this;
}


/**
 * Determines the maximum depth, in terms of tree levels, of the (sub) tree with this node as root
 *
 * @return The maximum number of tree levels of the subtree with this node as root
 */
int BTree::ParallelNode::maximumRecursionDepth() const {
        int maxChilds =
                std::max(
                    std::max(
                        std::max(
                                octNodes_[0]==nullptr ? 0 : octNodes_[0]->maximumRecursionDepth(),
                                octNodes_[1]==nullptr ? 0 : octNodes_[1]->maximumRecursionDepth()),
                        std::max(
                                octNodes_[2]==nullptr ? 0 : octNodes_[2]->maximumRecursionDepth(),
                                octNodes_[3]==nullptr ? 0 : octNodes_[3]->maximumRecursionDepth())),
                    std::max(
                        std::max(
                                octNodes_[4]==nullptr ? 0 : octNodes_[4]->maximumRecursionDepth(),
                                octNodes_[5]==nullptr ? 0 : octNodes_[5]->maximumRecursionDepth()),
                        std::max(
                                octNodes_[6]==nullptr ? 0 : octNodes_[6]->maximumRecursionDepth(),
                                octNodes_[7]==nullptr ? 0 : octNodes_[7]->maximumRecursionDepth()))
                );
        return 1 + maxChilds;
}


/**
 * Count nodes on individual tree levels in the subtree with this node as root.
 *
 * @param level The level of this node is on the global tree
 * @param numOfNodesOnLevels a vector, one element per tree level, to count the
 * the nodes number into
 */
void BTree::ParallelNode::countNodesOnLevel(int level, std::vector<int> &numOfNodesOnLevels) const {
            if(this->isRoot())
                numOfNodesOnLevels[0]++;
            for(auto octNode : octNodes_)
            {
                if(octNode!=nullptr)
                {
                    octNode->countNodesOnLevel(level+1, numOfNodesOnLevels);
                    numOfNodesOnLevels[level]++;
                }
            }
}


/**
 * Serializes the subtree with this node as root into a serialized vector of pointers to the nodes. The nodes are
 * sorted according to their level.
 * The current insert position on the individual tree levels are given in a second vector.
 * The given serializedNodes vector has to be sufficiently large to hold all nodes of the sub tree and the
 * initial insertPositions have to be at the correct positions so that nodes of the individual tree levels fit into
 * the serialized vector without overlap
 *
 * @param serializedNodes A vector of node pointers in which the subtree is serialized into
 * @param treeLevel The current tree level (starting from the root)
 * @param insertPositions A vector with the current insert positions in the serialized vector for the individual
 * tree levels
 */
void BTree::ParallelNode::serializeIntoVector(std::vector<BTree::ParallelNode*>& serializedNodes, int treeLevel,
                                              std::vector<int>& insertPositions)
{
    // Also updates the normalized edge length of the node, which is used for the
    // decision in the electric field calculation if the node is resolved in its subnodes
    updateNormalizedEdgeLength_();

    serializedNodes.at(insertPositions[treeLevel])=this;
    insertPositions[treeLevel]++;
    for(int i=0; i<8; ++i)
    {
        if(octNodes_[i]!=nullptr)
            octNodes_[i]->serializeIntoVector(serializedNodes, treeLevel+1, insertPositions);
    }
}


/**
 * Updates the edge length member according to the actual dimensions of the node
 */
void BTree::ParallelNode::updateNormalizedEdgeLength_(){
    edgeLengthSquaredNormalized_= ((max_.x()-min_.x()) * (max_.x()-min_.x())) / (BTree::AbstractNode::theta*BTree::AbstractNode::theta);
}
