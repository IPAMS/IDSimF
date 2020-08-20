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
 RS_ReactiveParticle.hpp

 A reactive version of a simulated particle

 ****************************/


#ifndef RS_ReactiveParticle_hpp
#define RS_ReactiveParticle_hpp

#include "BTree_particle.hpp"
#include "RS_Substance.hpp"

namespace RS{class ReactiveParticle;}
std::ostream& operator<<(std::ostream& os, const RS::ReactiveParticle& particle);

namespace RS {
    /**
     * A reactive particle: A particle (most probably ionic) which is able to react in a chemical reaction
     */
    class ReactiveParticle : public BTree::Particle {

    public:
        explicit ReactiveParticle(Substance* species);
        ReactiveParticle(Substance* species,Core::Vector location);
        ReactiveParticle(Substance* species,Core::Vector location, double charge);

        void setSpecies(Substance* species);
        Substance* getSpecies();

        friend std::ostream& ::operator<<(std::ostream& os, const RS::ReactiveParticle& particle);

    private:
        Substance *species_; ///< a link to a chemical substance this particle is made of
        void updateParticleParametersFromSpecies_();
    };
}

typedef std::unique_ptr<RS::ReactiveParticle> uniqueReactivePartPtr;

#endif //RS_ReactiveParticle_hpp
