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
 BTree_simpleVTKwriter.hpp

 File writer to write particle cloud and Barnes Hut tree data to legacy VTK files

 ****************************/

#ifndef BTree_simpleVTKwriter_hpp
#define BTree_simpleVTKwriter_hpp

#include <fstream>

//forward declare own classes:
namespace BTree{
    class Tree;
    class Node;
}

namespace ParticleSimulation {
    /**
     * File writer to write particle cloud and Core data to legacy VTK files
     */
    class SimpleVTKwriter {

    public:
        SimpleVTKwriter(std::string basefilename);
        ~SimpleVTKwriter();
        void write(BTree::Tree& tree,bool writeTree);
        
    private:
        std::ofstream* particleFile_; ///< filestream to write particle data to
        std::ofstream* treeFile_; ///< filestream to write tree data to

        static void writeIonPosition(std::ostream& filestream, const BTree::Node* rootnode);
        static void writeIonActive(std::ostream& filestream, const BTree::Node* rootnode);
        static void writeTreeBlockCornerCoordinates(std::ostream& filestream, const BTree::Node* rootnode);
        static void writeTreeBlockPolygon(std::ostream& filestream, int startIndex);
    };
}

#endif /* BTree_simpleVTKwriter_hpp */
