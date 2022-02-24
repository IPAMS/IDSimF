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

#include "FileIO_simpleVTKwriter.hpp"
#include "Core_particle.hpp"
#include "BTree_tree.hpp"
#include "BTree_node.hpp"
//#include <iostream>

/**
 * Constructor: Creates a VTK filewriter to write to result files.
 *
 * Two files are created:
 * basefilename_particles.vtk - which the particle trajectory data is written to and
 * basefilename_tree.vtk - which the tree / tree node data is written to
 *
 * @param basefilename basename for the result files to create
 */
FileIO::SimpleVTKwriter::SimpleVTKwriter(std::string basefilename){
    particleFile_ = new std::ofstream();
    treeFile_ = new std::ofstream();
    
    particleFile_->open(basefilename+"_particles.vtk");
    treeFile_->open(basefilename+"_tree.vtk");
}

/**
 * Destructor
 */
FileIO::SimpleVTKwriter::~SimpleVTKwriter(){
    particleFile_->flush();
    particleFile_->close();
    delete (particleFile_);
    
    treeFile_->flush();
    treeFile_->close();
    delete (treeFile_);
}

/**
 * Writes the current state of a Core to the result files
 * @param tree the Core to write to the file
 * @param writeTree if true: also tree node information is written to the tree file, otherwise only particle information
 * is written
 */
void FileIO::SimpleVTKwriter::write(BTree::Tree& tree,bool writeTree){
    //the base mechanism here is:
    //Core is able to export internal information through a given function to a file stream via the "writeToStream"
    //method
    //We provide the static data transform functions and pass them to the tree via function pointers and the
    //tree applies them to itself

    BTree::Node* root = tree.getRoot();
    std::size_t nIons = root->getNumberOfParticles();

    std::string header = "# vtk DataFile Version 2.0\n"
                         "Core Test\n"
                         "ASCII\n"
                         "DATASET POLYDATA\n"
                         "POINTS ";
    *particleFile_<<header<<nIons<<" float\n";
    root->writeToStream(*particleFile_, &SimpleVTKwriter::writeIonPosition); //pass the particle position function to the tree
    *particleFile_<<"POINT_DATA "<<nIons<<"\n";
    *particleFile_<<"SCALARS active integer 1\nLOOKUP_TABLE default\n";
    root->writeToStream(*particleFile_,&SimpleVTKwriter::writeIonActive); //pass the "particle active" function to the tree
    particleFile_->flush();

    if (writeTree){
        int nNodes = root->getNumberOfNodes();
        *treeFile_<<header<<nNodes*8<<" float\n";

        //pass node geometry write function to the tree:
        root->writeToStream(*treeFile_,&SimpleVTKwriter::writeTreeBlockCornerCoordinates);
        *treeFile_<<"POLYGONS "<<nNodes*6<<" "<<nNodes*30<<"\n";
        for (int i=0; i<nNodes;i++){
            writeTreeBlockPolygon(*treeFile_, i*8);
        }
        treeFile_->flush();
    }
}

/**
 * A function to export ion position information from a tree
 * @param filestream the filestream to write to
 * @param rootnode root node of the tree
 */
void FileIO::SimpleVTKwriter::writeIonPosition(std::ostream& filestream, const BTree::Node* rootnode){
    if (rootnode->getNumberOfParticles() == 1){
        BTree::TreeParticle* ion = rootnode->getParticle();
        filestream<<ion->wrappedParticle->getLocation()<<"\n";
    }
}

/**
 * A function to export information if ions are active from a tree
 * @param filestream the filestream to write to
 * @param rootnode root node of the tree
 */
void FileIO::SimpleVTKwriter::writeIonActive(std::ostream& filestream, const BTree::Node* rootnode){
    if (rootnode->getNumberOfParticles() == 1){
        BTree::TreeParticle* ion = rootnode->getParticle();
        if (ion->wrappedParticle->isActive() == true){
            filestream<<"1\n";
        }
        else{
            filestream<<"0\n";
        }
    }
}

/**
 * A function to export the corner coordinates of the tree nodes from a tree
 * @param filestream the filestream to write to
 * @param rootnode root node of the tree
 * */
void FileIO::SimpleVTKwriter::writeTreeBlockCornerCoordinates(std::ostream& filestream, const BTree::Node* rootnode){
    Core::Vector min= rootnode->getMin();
    Core::Vector max= rootnode->getMax();
    filestream<<min.x()<<" "<<min.y()<<" "<<min.z()<<"\n";
    filestream<<max.x()<<" "<<min.y()<<" "<<min.z()<<"\n";
    filestream<<max.x()<<" "<<max.y()<<" "<<min.z()<<"\n";
    filestream<<min.x()<<" "<<max.y()<<" "<<min.z()<<"\n";
    filestream<<min.x()<<" "<<min.y()<<" "<<max.z()<<"\n";
    filestream<<max.x()<<" "<<min.y()<<" "<<max.z()<<"\n";
    filestream<<max.x()<<" "<<max.y()<<" "<<max.z()<<"\n";
    filestream<<min.x()<<" "<<max.y()<<" "<<max.z()<<"\n";
}

/**
 * Writes the polygon node links for a tree node block to a filestream
 * @param filestream the filestream to write to
 * @param startIndex the index in the coodinate list of the first coordinate of the block to write links for
 */
void FileIO::SimpleVTKwriter::writeTreeBlockPolygon(std::ostream& filestream, int startIndex){
    //filestream<<startIndex<<"\n";
    filestream<<"4 "<<0+startIndex<<" "<<0+startIndex<<" "<<0+startIndex<<" "<<0+startIndex<<"\n";
    filestream<<"4 "<<4+startIndex<<" "<<5+startIndex<<" "<<6+startIndex<<" "<<7+startIndex<<"\n";
    filestream<<"4 "<<0+startIndex<<" "<<1+startIndex<<" "<<5+startIndex<<" "<<4+startIndex<<"\n";
    filestream<<"4 "<<2+startIndex<<" "<<3+startIndex<<" "<<7+startIndex<<" "<<6+startIndex<<"\n";
    filestream<<"4 "<<0+startIndex<<" "<<4+startIndex<<" "<<7+startIndex<<" "<<3+startIndex<<"\n";
    filestream<<"4 "<<1+startIndex<<" "<<2+startIndex<<" "<<6+startIndex<<" "<<5+startIndex<<"\n";
}
