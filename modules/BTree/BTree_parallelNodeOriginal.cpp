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

#include "BTree_parallelNodeOriginal.hpp"
#include <math.h>
#include <sstream>
#include <iostream>


//double Core::Node::theta = 0.0;

/* Class Constructor.
 param min the "lower" position of the node (one corner of the spatial block)
 param max the "upper" position of the node (the other corner of the spatial block)
 param parent the parenting node in the Core
 */
BTree::ParallelNodeOriginal::ParallelNodeOriginal(Core::Vector min, Core::Vector max, BTree::ParallelNodeOriginal* parent)
:
AbstractNode(min,max),
parent_(parent)
{
    for (int i=0; i<8; i++){
        this->octNodes_[i] = nullptr;
    }
    nNodes_++;
}



// 1. copy constructor
/*person(const person& that)
{
    name = new char[strlen(that.name) + 1];
    strcpy(name, that.name);
    age = that.age;
}*/

BTree::ParallelNodeOriginal::ParallelNodeOriginal(const BTree::ParallelNodeOriginal& that):
AbstractNode(that),
parent_(that.parent_)
{
    for (int i=0; i<8; i++){
        this->octNodes_[i] = that.octNodes_[i];
    }
    nNodes_++;
}

BTree::ParallelNodeOriginal& BTree::ParallelNodeOriginal::operator=(const BTree::ParallelNodeOriginal& that){

    max_=that.max_;
    min_=that.min_;
    center_=that.center_;
    centerOfCharge_=that.centerOfCharge_;
    charge_=that.charge_;
    parent_=that.parent_;
    particle_=that.particle_;
    numP_=that.numP_;

    for (int i=0; i<8; i++){
        this->octNodes_[i] = that.octNodes_[i];
    }
    nNodes_++;

    return *this;
}

/* Class Destructor. 
 destroys the node
 */
BTree::ParallelNodeOriginal::~ParallelNodeOriginal(){
    //FIXME std::cout<<"node destroyed "<<this<<std::endl;
    for (int i=0; i<8; i++){
        if(octNodes_[i] != nullptr){
            delete(octNodes_[i]);
            octNodes_[i] = nullptr;
        }
    }
    nNodes_--;
}

// Static methods:
/**
 * Inits a node as root of a tree
 */
void BTree::ParallelNodeOriginal::initAsRoot(){
    //FIXME: Check if the current parent has to be set to nullptr
    //Core::Node::nNodes_=1; //at least the root node has to exist
    this->parent_ = nullptr;
}

/**
 * Returns true if the node is a root (has no parent)
 */
bool BTree::ParallelNodeOriginal::isRoot() const{
    if (parent_ == nullptr){
        return true;
    }
    else {
        return false;
    }
}

/* todo: comment */
BTree::ParallelNodeOriginal::Octant BTree::ParallelNodeOriginal::getOctant(Core::Vector location){
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


/**
 * Inserts a particle into the node. If the node is already occupied, the node is subdivided into sub-nodes and
 * the particle is inserted in the according sub-nodes.
 *
 * @param particle the particle to insert into the node.
 */
void BTree::ParallelNodeOriginal::insertParticle(BTree::Particle* particle){

    if (numP_ >1){
        Octant oct = this->getOctant(particle->getLocation());
        if (this->octNodes_[oct] == nullptr){
            this->octNodes_[oct] = this->createOctNode(oct);
        }
        octNodes_[oct]->insertParticle(particle);
    }
    else if(numP_ == 1){
        BTree::Particle* p2 = particle_;
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


void BTree::ParallelNodeOriginal::removeMyselfFromTree(){

    //iterate through parent nodes:
    BTree::ParallelNodeOriginal* parentNode = parent_;
    BTree::ParallelNodeOriginal* leafNode = nullptr;
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

void BTree::ParallelNodeOriginal::updateSelf(){
    centerOfCharge_ = particle_->getLocation();
    charge_ = particle_->getCharge();
}

void BTree::ParallelNodeOriginal::updateParents(){
    BTree::ParallelNodeOriginal* parentNode = this->parent_;
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


/* todo: comment */
BTree::ParallelNodeOriginal* BTree::ParallelNodeOriginal::createOctNode(BTree::ParallelNodeOriginal::Octant oct){

    BTree::ParallelNodeOriginal* result = nullptr;
    if (oct == SWB){
        result = new ParallelNodeOriginal(min_,center_,this);
    }
    else if(oct == NWB){
        result = new ParallelNodeOriginal(
                Core::Vector(min_.x(),center_.y(),min_.z()),
                Core::Vector(center_.x(),max_.y(),center_.z()),
                this);
    }
    else if(oct == NEB){
        result = new ParallelNodeOriginal(
                Core::Vector(center_.x(),center_.y(),min_.z()),
                Core::Vector(max_.x(),max_.y(),center_.z()),
                this);
    }
    else if(oct == SEB){
        result = new ParallelNodeOriginal(
                Core::Vector(center_.x(),min_.y(),min_.z()),
                Core::Vector(max_.x(),center_.y(),center_.z()),
                this);
    }
    else if(oct == SWT){
        result = new ParallelNodeOriginal(
                Core::Vector(min_.x(),min_.y(),center_.z()),
                Core::Vector(center_.x(),center_.y(),max_.z()),
                this);
    }
    else if(oct == NWT){
        result = new ParallelNodeOriginal(
                Core::Vector(min_.x(),center_.y(),center_.z()),
                Core::Vector(center_.x(),max_.y(),max_.z()),
                this);
    }
    else if(oct == NET){
        //result = new Node(center_,max_,this);
        result = new ParallelNodeOriginal(
                Core::Vector(center_.x(),center_.y(),center_.z()),
                Core::Vector(max_.x(),max_.y(),max_.z()),
                this);
    }
    else if(oct == SET){
        result = new ParallelNodeOriginal(
                Core::Vector(center_.x(),min_.y(),center_.z()),
                Core::Vector(max_.x(),center_.y(),max_.z()),
                this);
    }

    return result;
}



void BTree::ParallelNodeOriginal::setDurchmesser(){
    // Gleichung wurde umgestellt. In setDurchmesser wird jetzt der Durchmesser zum Quadrat berechnet und druch Theta² geteilt
    durch_=((max_.x()-min_.x())*(max_.x()-min_.x()))/(BTree::AbstractNode::theta*BTree::AbstractNode::theta);
}

// Public member methods:

int BTree::ParallelNodeOriginal::max_rek_tiefe(BTree::ParallelNodeOriginal* nod)
{
    if(nod==nullptr)
        return 0;
    else{
        int maxChilds = +std::max(std::max(std::max(max_rek_tiefe(nod->octNodes_[0]),max_rek_tiefe(nod->octNodes_[1])),std::max(max_rek_tiefe(nod->octNodes_[2]),max_rek_tiefe(nod->octNodes_[3]))),std::max(std::max(max_rek_tiefe(nod->octNodes_[4]),max_rek_tiefe(nod->octNodes_[5])),std::max(max_rek_tiefe(nod->octNodes_[6]),max_rek_tiefe(nod->octNodes_[7]))));
        //FIXME std::cout<<" "<<maxChilds<<std::endl;
        return 1 + maxChilds;
    }

}

void BTree::ParallelNodeOriginal::bestimme_Anfang(std::vector<int> &begin, std::vector<int> &zaehler, int n)
{
    begin[0]=0;
    for(int i=1;i<n;i++)
    {
        begin.at(i)=begin.at(i-1)+zaehler.at(i-1);
    }
}

void BTree::ParallelNodeOriginal::anz_Knoten_auf_der_Ebene(BTree::ParallelNodeOriginal* node, int ebene, std::vector<int> &zaehler, int n)
{
            if(node->isRoot())
                zaehler[0]++;
            for(int i=0;i<8;i++)
            {
                if(node->octNodes_[i]!=nullptr)
                {
                    anz_Knoten_auf_der_Ebene(node->octNodes_[i],ebene+1,zaehler,n);
                    zaehler[ebene]++;
                }
            }
}

void BTree::ParallelNodeOriginal::insertNodes(BTree::ParallelNodeOriginal* nod, std::vector<BTree::ParallelNodeOriginal*> &Nodes, int ebene, std::vector<int> &zaehler, std::vector<int> &begin, int n)
{
    // berechnet den Durchmesser für jeden Knoten
    nod->setDurchmesser();
    Nodes.at(begin[ebene-1])=nod;
    begin[ebene-1]++;
    for(int i=0;i<8;i++)
    {
        if(nod->octNodes_[i]!=nullptr)
            insertNodes(nod->octNodes_[i],Nodes,ebene+1,zaehler,begin,n);
    }
}

void BTree::ParallelNodeOriginal::update(std::vector<BTree::ParallelNodeOriginal*> &Nodes, std::vector<int> &begin, int n, int summe)
{
    //FIXME: std::cout << "n1: "<<n<<std::endl;
    int m=begin[n-1];
    int j;
    #pragma omp parallel for
    for(int k=summe-1;k>=m;k--)
    {
        if(Nodes.at(k) != nullptr){
                if(Nodes.at(k)->numP_==1)
                {
                    Nodes[k]->centerOfCharge_ = Nodes[k]->particle_->getLocation();
                    Nodes[k]->charge_ = Nodes[k]->particle_->getCharge();
                }
                else
                {
                    
                    Nodes[k]->charge_ = 0.0;
                    Nodes[k]->centerOfCharge_ = Core::Vector(0.0,0.0,0.0);
                        for (int h=0; h<8; h++){
                            if(Nodes[k]->octNodes_[h] != nullptr){
                                Nodes[k]->charge_ += Nodes[k]->octNodes_[h]->getCharge();
                                Nodes[k]->centerOfCharge_ = Nodes[k]->centerOfCharge_+
                                (Nodes[k]->octNodes_[h]->centerOfCharge_ *
                                Nodes[k]->octNodes_[h]->charge_);
                            }
                        }
            
                   if (Nodes[k]->charge_ != 0.0){
                    Nodes[k]->centerOfCharge_ = Nodes[k]->centerOfCharge_ / Nodes[k]->charge_;
                  }
                  else {
                    Nodes[k]->centerOfCharge_ = Nodes[k]->center_;
                  }
                }
        }
    }

    //FIXME: std::cout << "n2: "<<n<<std::endl;
    for(int i=n;i>=2;i--)
    {
        int k=begin.at(i-1);
        int l=begin.at(i-2);
        #pragma omp parallel for
        for(j=k-1;j>=l;j--)
        {
            if(Nodes.at(j) != nullptr){
                if(Nodes.at(j)->numP_==1)
                {
                    Nodes[j]->centerOfCharge_ = Nodes[j]->particle_->getLocation();
                    Nodes[j]->charge_ = Nodes[j]->particle_->getCharge();
                }
                else
                {
                    Nodes[j]->charge_ = 0.0;
                    Nodes[j]->centerOfCharge_ = Core::Vector(0.0,0.0,0.0);
                        for (int m=0; m<8; m++){
                            if(Nodes[j]->octNodes_[m] != nullptr){
                                Nodes[j]->charge_ += Nodes[j]->octNodes_[m]->getCharge();
                                Nodes[j]->centerOfCharge_ = Nodes[j]->centerOfCharge_+
                                (Nodes[j]->octNodes_[m]->centerOfCharge_ *
                                Nodes[j]->octNodes_[m]->charge_);
                            }
                        }
            
                    if (Nodes[j]->charge_ != 0.0){
                        Nodes[j]->centerOfCharge_ = Nodes[j]->centerOfCharge_ / Nodes[j]->charge_;
                    }
                    else {
                        Nodes[j]->centerOfCharge_ = Nodes[j]->center_;
                    }
                }
            }
        }
    }
}

Core::Vector BTree::ParallelNodeOriginal::computeEFieldFromTree(BTree::Particle& targetP, std::vector<BTree::ParallelNodeOriginal*> &MyNodes){
    //Core::Vector efield;

    //todo: Check if i caluclate the field to myself...
    if (numP_ == 1){
        return(calculateElectricField(
                targetP.getLocation(),
                particle_->getLocation(),
                particle_->getCharge()));
    }
    else{
        // magnitude2 berechnet den Abstand der beiden Werte zum Quadrat. Es wird nicht mehr die Wurzel gezogen
        double r = (targetP.getLocation()-centerOfCharge_).magnitudeSquared();
        // durch_ wird in der insert Funktion gesetzt und wird hier nur aufgerufen
        // durch_ wird anders berechnet als vorher. Theta wurde bei der Rechnung hinzugefügt, sodass die Ungleichung umgestellt wurde
        if(durch_ < r){
            return(calculateElectricField(
                    targetP.getLocation(),
                    centerOfCharge_,
                    charge_));

        }
        else{
            Core::Vector efield= Core::Vector(0.0,0.0,0.0);
            Core::Vector loc=targetP.getLocation();
            // Zähler für die aktuelle Position im Feld
            int j=0;
            // Zähler für die nächste freie Position im Feld
            int k=1;

            // fügt den aktuellen Knoten in die aktuelle Position ein
            MyNodes.at(j) = this;
            while(j<k)
            {
                // Berechnet die Kraft explizit für jeden Knoten der in das Feld hinzugefügt wurde.
                // Hier wird der erste Fall umgesetzt. Sollte der aktuelle Knoten nur ein Partikel beinhalten
                if(MyNodes[j]->numP_==1)
                    efield=efield+calculateElectricField(loc, MyNodes[j]->particle_->getLocation(),
                            MyNodes[j]->particle_->getCharge());
                // Ansonsten wird die Entfernung zum Masseschwerpunkt bestimmt
                else
                {
                    double r2= (loc-MyNodes[j]->centerOfCharge_).magnitudeSquared();
                    if(MyNodes[j]->durch_<r2)
                        efield=efield+calculateElectricField(loc, MyNodes[j]->centerOfCharge_, MyNodes[j]->charge_);
                    // Sollte die Entfernung größer als r2 sein, dann werden die noch zu betrachtenden Knoten dem Feld hinzugefügt
                    else
                    {
                        for(int i=0; i<8; i++)
                        {
                            if(MyNodes[j]->octNodes_[i]!=nullptr)
                            {
                              // fügt den nachfolgenden Knoten an die nächste freie Stelle
                              MyNodes[k]=MyNodes[j]->octNodes_[i];
                              // erhöht den Zähler für die nächste freie Stelle um eins
                              k++;
                            }
                        }
                    }
                }
                // j wird am Ende der Schleife immer um eins erhöht.
                // Sollte in der for Schleife kein neuer Knoten hinzugekommen sein, dann ist j=k und die Schleife wird abgebrochen
                j++;
            }
            return(efield);
        }
    }
}

Core::Vector BTree::ParallelNodeOriginal::computeElectricForceFromTree(BTree::Particle&){
    std::stringstream ss;
    ss << "Method not implemented: Core::ParallelNodeOriginal::computeElectricFieldFromTree";
    throw (std::runtime_error(ss.str()));
}

//Todo: comment method
void BTree::ParallelNodeOriginal::computeChargeDistributionRecursive(){
    if (numP_ == 1){
        //this->updateSelf();
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
std::string BTree::ParallelNodeOriginal::toString() const{
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
void BTree::ParallelNodeOriginal::printTree(int level) const{

    for (int i=0; i<8; i++){
        if(octNodes_[i] != nullptr){
            octNodes_[i]->printTree(level+1);
        }
    }

    for (int i=0; i<level; i++){
        std::cout<<"--|";
    }
    std::cout<<"this:"<<this<<" min "<<this->min_<<" max "<<this->max_<<" part "<<this->numP_<<" charge "<<this->charge_;
    //std::cout<<" min "<<this->min_<<" max "<<this->max_<<" part "<<this->numP_<<" charge "<<this->charge_;

    if (this->particle_ != nullptr){
        std::cout<<" pl:"<<this->particle_<<" loc:"<<this->particle_->getLocation();
    }
    std::cout<<std::endl;
}

void BTree::ParallelNodeOriginal::writeToStream(std::ostream& filestream,
                                                void (*writeFct)(std::ostream& filestream, const BTree::ParallelNodeOriginal* node) ) const{
    
    writeFct(filestream,this);
    
    for (int i=0; i<8; i++){
        if(octNodes_[i] != nullptr){
            octNodes_[i]->writeToStream(filestream,writeFct);
            
        }
    }
}


