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

 ------------
 BTree_genericBaseNode.hpp

 Templated base class for nodes in Barnes-Hut tree, which bundles generalizable
 methods / algorithms for concrete tree node types NodType
 ****************************/

//FIXME: put template implementation to cpp file and use explicit instantiation
#ifndef BTree_generic_base_node_hpp
#define BTree_generic_base_node_hpp


#include "Core_vector.hpp"
#include "BTree_abstractNode.hpp"
#include "BTree_particle.hpp"
#include <sstream>
#include <iostream>

namespace BTree{

    /**
     * Templated base class for nodes in Barnes-Hut tree, which bundles generalizable
     * methods / algorithms for concrete tree node types NodType.
     *
     * This class uses the Curiously recurring template pattern (CRTP) / F-bound idiom to generalize
     * functionality for concrete tree nodes.
     *
     * @tparam NodType the concrete node type to generate generalizable methods for
     */
    template<class NodType>
    class GenericBaseNode: public AbstractNode {

    public:

        //Constructors:
        GenericBaseNode(Core::Vector min, Core::Vector max, NodType* parent);
        GenericBaseNode(const GenericBaseNode& that) = delete;
        GenericBaseNode& operator=(const GenericBaseNode& that) = delete;

        // Destructor:
        ~GenericBaseNode() override;

        // Generic node / tree management methods:
        void initAsRoot() override;
        bool isRoot() const override;
        NodType* createOctNode(Octant oct);
        NodType** getOctants();
        void removeMyselfFromTree() override;
        void updateSelf() override;
        void updateParents() override;
        void insertParticle(Particle* particle) override;
        void computeChargeDistributionRecursive() override;

        // Diagnostic methods:
        [[nodiscard]] std::string toString() const override;
        virtual void printTree(int level) const;
        void writeToStream(std::ostream& filestream,void (*writeFct)(std::ostream& filestream, const NodType* node)) const;
        virtual void testNodeIntegrity(int level);


    protected:
        NodType* octNodes_[8]; ///< Subnodes of this node
        NodType* parent_; ///< Parent node of this node
    };


    /**
     * Basic constructor
     * @param min The lower corner of the node
     * @param max The upper corner of the node
     * @param parent A node which is the parent node in the tree structure of this node
     */
    template<class NodType>
    GenericBaseNode<NodType>::GenericBaseNode(Core::Vector min, Core::Vector max, NodType* parent):
        AbstractNode(min,max),
        parent_(parent)
    {
        for (int i=0; i<8; i++){
            this->octNodes_[i] = nullptr;
        }
        nNodes_++;
    }

    /**
     * Destructor
     */
    template<class NodType>
    GenericBaseNode<NodType>::~GenericBaseNode<NodType>(){
        //destroy all sub-nodes:
        for (int i=0; i<8; i++){
            if(octNodes_[i] != nullptr){
                delete(octNodes_[i]);
                octNodes_[i] = nullptr;
            }
        }
        nNodes_--;
    }


    /**
     * Inits a node as root of a tree
     */
    template<class NodType>
    void GenericBaseNode<NodType>::initAsRoot(){
        this->parent_ = nullptr;
    }

    /**
     * Returns true if the node is a root (has no parent)
     */
    template<class NodType>
    bool GenericBaseNode<NodType>::isRoot() const{
        if (parent_ == nullptr){
            return true;
        }
        else {
            return false;
        }
    }

    /**
     * Creates a new sub-octant node as a child of this node in one of the sub octants of this node.
     *
     * @param oct the location / sub-octant where the new node should be created
     * @return the newly created sub-octant node
     */
    template<class NodType>
    NodType* GenericBaseNode<NodType>::createOctNode(Octant oct) {

        NodType* parent = static_cast<NodType*>(this);

        NodType* result = nullptr;
        if (oct == SWB){
            result = new NodType(min_,center_,parent);
        }
        else if(oct == NWB){
            result = new NodType(
                    Core::Vector(min_.x(),center_.y(),min_.z()),
                    Core::Vector(center_.x(),max_.y(),center_.z()),
                    parent);
        }
        else if(oct == NEB){
            result = new NodType(
                    Core::Vector(center_.x(),center_.y(),min_.z()),
                    Core::Vector(max_.x(),max_.y(),center_.z()),
                    parent);
        }
        else if(oct == SEB){
            result = new NodType(
                    Core::Vector(center_.x(),min_.y(),min_.z()),
                    Core::Vector(max_.x(),center_.y(),center_.z()),
                    parent);
        }
        else if(oct == SWT){
            result = new NodType(
                    Core::Vector(min_.x(),min_.y(),center_.z()),
                    Core::Vector(center_.x(),center_.y(),max_.z()),
                    parent);
        }
        else if(oct == NWT){
            result = new NodType(
                    Core::Vector(min_.x(),center_.y(),center_.z()),
                    Core::Vector(center_.x(),max_.y(),max_.z()),
                    parent);
        }
        else if(oct == NET){
            //result = new Node(center_,max_,this);
            result = new NodType(
                    Core::Vector(center_.x(),center_.y(),center_.z()),
                    Core::Vector(max_.x(),max_.y(),max_.z()),
                    parent);
        }
        else if(oct == SET){
            result = new NodType(
                    Core::Vector(center_.x(),min_.y(),center_.z()),
                    Core::Vector(max_.x(),center_.y(),max_.z()),
                    parent);
        }

        return result;
    }


    /**
     * Get links to the child octant nodes
     * @return c style array with pointers to the child octant nodes
     */
    template<class NodType>
    NodType** GenericBaseNode<NodType>::getOctants(){
        return octNodes_;
    }


    /**
     * Removes this node from the tree structure
     */
    template<class NodType>
    void GenericBaseNode<NodType>::removeMyselfFromTree() {
        //iterate through parent nodes:
        NodType* parentNode = parent_;
        NodType* leafNode = nullptr;
        int lnodes=0;

        //remove this node from the octant nodes of the parent node:
        if (parentNode != nullptr){
            Octant oct = parentNode->getOctant(this->particle_->getLocation());

            //delete(parentNode->octNodes_[oct]);
            parentNode->octNodes_[oct] = nullptr;
        }
        else {
            //special condition: we are in the root: empty root node
            this->numP_ = this->numP_-1;
            this->particle_ = nullptr;
            this->charge_ = 0.0;
            this->centerOfCharge_ = Core::Vector(0,0,0);
        }

        //update the parent nodes up to the root:
        while (parentNode != nullptr) {
            parentNode->numP_ = parentNode->numP_ -1;

            //if only one particle is left in a parenting node, this node has no
            //childs anymore: it becomes the new leaf
            if (parentNode->numP_ == 1){
                for (int i=0; i<8; i++){
                    if(parentNode->octNodes_[i] != nullptr){

                        parentNode->particle_ = parentNode->octNodes_[i]->particle_;
                        parentNode->charge_ = parentNode->octNodes_[i]->charge_;
                        parentNode->particle_->setHostNode(parentNode);
                        leafNode = parentNode;
                        lnodes++;

                        //delete the sibling from the tree:
                        delete(parentNode->octNodes_[i]);
                        parentNode->octNodes_[i] = nullptr;
                        break;
                    }
                }
            }
            parentNode = parentNode->parent_;
        }

        if (leafNode == nullptr){
            this->updateParents();
        }
        else {
            leafNode->updateSelf();
            leafNode->updateParents();
        }
    }


    /**
     * Updates the state of this node
     */
    template<class NodType>
    void GenericBaseNode<NodType>::updateSelf(){
        centerOfCharge_ = particle_->getLocation();
        charge_ = particle_->getCharge();
    }


    /**
     * Updates the internal state of the parent nodes of this node up to the tree root
     */
    template<class NodType>
    void GenericBaseNode<NodType>::updateParents() {
        NodType* parentNode = this->parent_;
        while (parentNode != nullptr) {
            parentNode->charge_ = 0.0;
            parentNode->centerOfCharge_ = Core::Vector(0.0,0.0,0.0);
            for (int i=0; i<8; i++){
                if(parentNode->octNodes_[i] != nullptr){
                    parentNode->charge_ += parentNode->octNodes_[i]->getCharge();
                    parentNode->centerOfCharge_ = parentNode->centerOfCharge_+
                                                  (parentNode->octNodes_[i]->centerOfCharge_ *
                                                   parentNode->octNodes_[i]->charge_);
                }
            }
            if (parentNode->charge_ != 0.0){
                parentNode->centerOfCharge_ = parentNode->centerOfCharge_ / parentNode->charge_;
            }
            else {
                parentNode->centerOfCharge_ = parentNode->center_;
            }
            parentNode = parentNode->parent_;
        }
    }


    /**
     * Inserts a particle into the node. If the node is already occupied, the node is subdivided into sub-nodes and
     * the particle is inserted in the according sub-nodes.
     *
     * @param particle the particle to insert into the node.
     */
    template<class NodType>
    void GenericBaseNode<NodType>::insertParticle(Particle* particle){

        if (numP_ >1){
            Octant oct = this->getOctant(particle->getLocation());
            if (this->octNodes_[oct] == nullptr){
                this->octNodes_[oct] = this->createOctNode(oct);
            }
            octNodes_[oct]->insertParticle(particle);
        }
        else if(numP_ == 1){
            Particle* p2 = particle_;
            if(p2->getLocation() != particle->getLocation()){
                //There is already a particle in the node
                //relocate and subdivide
                Octant oct = this->getOctant(p2->getLocation());
                if (octNodes_[oct] == nullptr){
                    octNodes_[oct] = this->createOctNode(oct);
                }
                octNodes_[oct]->insertParticle(p2);


                particle_= nullptr;

                oct = this->getOctant(particle->getLocation());
                if (octNodes_[oct] == nullptr){
                    octNodes_[oct] = this->createOctNode(oct);
                }
                octNodes_[oct]->insertParticle(particle);
            }
            else {
                //if two particles with exactly the same position are existing: Throw exception
                std::stringstream ss;
                ss << "Tried to insert particle with exactly the same position: "<<particle->getLocation()<<std::endl;
                throw (std::logic_error(ss.str()));
            }

        }else{
            particle_ = particle;
            particle->setHostNode(this);
            this->updateSelf();
        }
        numP_++;
    }


    /**
     * Computes the charge distribution in this node and its subnodes and updates the node and its
     * subnodes with updates charge and center of charge values
     */
    template<class NodType>
    void GenericBaseNode<NodType>::computeChargeDistributionRecursive(){
        if (numP_ == 1){
            this->updateSelf();
            //std::cout << " compCharg 01 charge "<<charge_<<" location "<<centerOfCharge_<<std::endl;
        }
        else{
            charge_ = 0.0;
            centerOfCharge_ = Core::Vector(0,0,0);
            for (int i=0; i<8; i++){
                if(octNodes_[i] != nullptr){
                    octNodes_[i]->computeChargeDistributionRecursive();
                    charge_ += octNodes_[i]->getCharge();
                    centerOfCharge_ = centerOfCharge_+
                            (octNodes_[i]->getCenterOfCharge() *
                                    octNodes_[i]->getCharge());
                }
            }
            if (charge_ != 0.0){
                centerOfCharge_ = centerOfCharge_ / charge_;
            }
            else {
                centerOfCharge_ = center_;
            }
        }
    }


    /**
     * Returns a string which reflects the state of the node
     */
    template<class NodType>
    std::string GenericBaseNode<NodType>::toString() const{
        std::stringstream ss;
        ss<<"this: "<<this<<" min "<<this->min_<<" max "<<this->max_<<" part "<<this->numP_<<" charge "<<this->charge_/Core::ELEMENTARY_CHARGE<<std::endl;

        if (this->particle_ != nullptr){
            ss<<"part: "<<this->particle_<<std::endl;
        }

        ss<<" >> childs: ";
        for (int i=0; i<8; i++){
            if(octNodes_[i] != nullptr){
                ss<<"i: "<<i<<" ch:"<<octNodes_[i]<<" nPar:"<<octNodes_[i]->numP_<<"| ";
            }
        }
        ss<<std::endl;

        if (this->particle_ != nullptr){
            ss<<"particle location:"<<this->particle_->getLocation()<<std::endl;
        }

        return(ss.str());
    }


    /**
     * Prints the tree with this node as root to cout
     * @param level level of the current tree node with respect to the tree root (should be 0 when called manually)
     */
    template<class NodType>
    void GenericBaseNode<NodType>::printTree(int level) const{

        for (int i=0; i<8; i++){
            if(octNodes_[i] != nullptr){
                octNodes_[i]->printTree(level+1);
            }
        }

        for (int i=0; i<level; i++){
            std::cout<<"--|";
        }
        std::cout<<" min "<<this->min_<<" max "<<this->max_<<" part "<<this->numP_<<" charge "<<this->charge_;

        if (this->particle_ != nullptr){
            std::cout<<" pl:"<<this->particle_<<" loc:"<<this->particle_->getLocation();
        }
        std::cout<<std::endl;
    }


    /**
     * A method to export internal information to a filestream. The data transformation and formatting is done
     * by a passed write function. The write function is recursively applied to all subnodes of this node and thus
     * the whole subtree below this node.
     *
     * @param filestream a filestream to write to
     * @param writeFct a function, taking the filestream and a Core node, which extracts information from the node
     */
    template<class NodType>
    void GenericBaseNode<NodType>::writeToStream(std::ostream& filestream,
                                                 void (*writeFct)(std::ostream& filestream, const NodType* node) ) const{

        writeFct(filestream, static_cast<const NodType*>(this));

        for (int i=0; i<8; i++){
            if(octNodes_[i] != nullptr){
                octNodes_[i]->writeToStream(filestream,writeFct);

            }
        }
    }


    /**
     * Tests the node integrity of the subtree with the current node as root
     *
     * @param debug prints some additional information to cout if true
     * @param level level of the current tree node with respect to the tree root (should be 0 when called manually)
     *
     * @throws logic_error if an inconsistent node is found
     */
    template<class NodType>
    void GenericBaseNode<NodType>::testNodeIntegrity(int level){

        if (this->numP_ == 0){
            std::stringstream ss;
            ss<<" >>>>> Zero particle node found :"<<std::endl;
            ss<<this->toString();
            throw (std::logic_error(ss.str()));
        }

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

        for (int i=0; i<8; i++){
            if(octNodes_[i] != nullptr){
                if (octNodes_[i]->parent_ != this){
                    std::stringstream ss;
                    ss << " >>>>> inconsistent node found this-min:" << this->min_ << " max " << this->max_ << " nPart"
                       << this->numP_ << " " << " sub-parent-min:" << octNodes_[i]->parent_->min_ << " max "
                       << octNodes_[i]->parent_->max_ << std::endl;
                    ss << "this:" << this->min_.x() - this->max_.x() << " parent:"
                       << octNodes_[i]->parent_->min_.x() - octNodes_[i]->parent_->max_.x() << std::endl;
                    ss << this << " " << &(octNodes_[i]) << " " << &octNodes_[i] << " " << octNodes_[i] << " "
                       << octNodes_[i]->parent_ << " i:" << i << " | " << &(octNodes_[i]->parent_) << std::endl;
                    throw (std::logic_error(ss.str()));
                }
                octNodes_[i]->testNodeIntegrity(level+1);
            }
        }
    }
}

#endif //BTree_generic_base_node_hpp
