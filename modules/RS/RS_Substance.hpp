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
 RS_Substance.hpp
 
 This implements a chemical substance in the RS simulation

 ****************************/

#ifndef RS_Substance_hpp
#define RS_Substance_hpp

#include <iostream>
#include <string>

namespace RS{class Substance;}
std::ostream& operator<<(std::ostream& os, const RS::Substance& subst);

namespace RS {

    // individual exception for problems in ion cloud files
    class SubstanceException : public std::exception
    {
    public:
        explicit SubstanceException (std::string s) {this->message = s;}
        ~SubstanceException() override = default;
        const char * what() const throw() override {return message.c_str();}

    private:
        std::string message;
    };

    class Substance {

    public:
        enum substanceType {
            isotropic, ///< isotropic continuous substance with isotropic concentration
            discrete, ///< substance modeled as individual simulated molecules
            field ///< non isotropic continuous substance with non isotropic concentration given as field
        };
        
        Substance(std::string name, substanceType type);
        Substance(std::string name, std::string typeLabel);

        [[nodiscard]] std::string name() const;
        [[nodiscard]] substanceType type() const;

        [[nodiscard]] double charge() const;
        void charge(double newCharge);
        [[nodiscard]] double mass() const;
        void mass(double newMass);
        [[nodiscard]] double mobility() const;
        void mobility(double newMobility);
        [[nodiscard]] double collisionDiameter() const;
        void collisionDiameter(double newCollisionDiameter);
        [[nodiscard]] double staticConcentration() const;
        void staticConcentration(double newStaticConcentration);

        friend std::ostream& ::operator<<(std::ostream& os, const RS::Substance& subst);

    private:
        std::string name_; ///< The name of this substance
        double staticConcentration_ = 0.0; ///< the static background concentration (for isotropic substances)
        double mass_ = 0.0; ///< the mass of molecules of this substance
        double charge_ = 0.0; ///< the charge of molecules of this substance
        double mobility_ = 0.0; ///< the electrical ion mobility of molecules of this substance
        double collisionDiameter_ = 0.0; ///< the effective collision diameter of molecules of this substance
        substanceType type_; ///< the type of this substance

        explicit Substance(std::string name);
    };

    bool operator<(const Substance& a, const Substance& b);
}

#endif /* RS_Substance_hpp */
