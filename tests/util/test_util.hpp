#ifndef IONSIMULATION_CPP_TEST_UTIL_HPP
#define IONSIMULATION_CPP_TEST_UTIL_HPP

#include "Core_vector.hpp"
#include "PSim_boxStartZone.hpp"
#include <string>
#include <sstream>
#include <fstream>
#include "catch.hpp"

const std::string vectorsApproxEqual("Vectors approximately equal");

inline std::string vectorApproxCompare(const Core::Vector &lhs, const Core::Vector &rhs){
    bool equal = true;
    std::stringstream ss;


    if (lhs.x() != Approx(rhs.x())){
        equal = false;
        ss << "x coord unequal "<<lhs.x() <<" != "<<rhs.x()<<std::endl;
    }

    if (lhs.y() != Approx(rhs.y())){
        equal = false;
        ss << "y coord unequal "<<lhs.y() <<" != "<<rhs.y()<<std::endl;
    }

    if (lhs.z() != Approx(rhs.z())){
        equal = false;
        ss << "z coord unequal "<<lhs.z() <<" != "<<rhs.z()<<std::endl;
    }

    if (equal){
        return vectorsApproxEqual;
    }
    else{
        return ss.str();
    }
}

//FIXME: get rid of this method: Rewrite tests with start zone directly
inline std::vector<std::unique_ptr<BTree::Particle>> getRandomIonsInBox(int numIons, Core::Vector corner, Core::Vector boxSize){
    Core::Vector center = corner + boxSize/2.0;
    ParticleSimulation::BoxStartZone startZone(boxSize, center);

    return startZone.getRandomParticlesInStartZone(numIons, 1.0);
}

inline std::string readTextFile(std::string filename){
    std::ifstream inputFilestream(filename);
    std::string result(
            (std::istreambuf_iterator<char>(inputFilestream)),
            (std::istreambuf_iterator<char>()));

    return result;
}

#endif //IONSIMULATION_CPP_TEST_UTIL_HPP
