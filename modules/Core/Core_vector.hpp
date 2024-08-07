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

 Core_vector.hpp
 
 This implements a three dimensional vector

 ****************************/

#ifndef BTree_vector_hpp
#define BTree_vector_hpp

#include <fstream>

namespace Core {

    /**
     * Simple 3d vector
     */
    class Vector {

    public:
        // Constructors:
        Vector() = default;
        Vector(double x, double y, double z);
        explicit Vector(const double*);
        
        // Accessors:
        [[nodiscard]] double x() const;
        [[nodiscard]] double y() const;
        [[nodiscard]] double z() const;
        
        // Setters:
        void set(double x,double y,double z);
        void x(double x);
        void y(double y);
        void z(double z);
        
        //mathematical methods:
        [[nodiscard]] double magnitude() const;
        [[nodiscard]] double magnitudeSquared() const;
        [[nodiscard]] Vector crossProduct(const Vector &rhs) const;

        // overloaded operators:
        Vector& operator+=(const Vector &rhs);
        friend Vector operator+(const Vector &lhs, const Vector &rhs);
        friend Vector operator-(const Vector &lhs, const Vector &rhs);
        friend double operator*(const Vector &lhs, const Vector &rhs);
        friend Vector operator*(const Vector &lhs, double rhs);
        friend Vector operator/(const Vector &lhs, double rhs);
        friend Vector operator*(double lhs, const Vector &rhs);

        friend bool operator==(Vector const &lhs, Vector const &rhs);
        friend bool operator!=(const Vector &lhs, const Vector &rhs);

        // Public member methods:
        void printState() const;

    private:
        // the x,y,z components of the vector
        double x_= 0.0;
        double y_= 0.0;
        double z_= 0.0;
    };


    // overloaded operators:
    Vector operator+(const Vector &lhs, const Vector &rhs);
    Vector operator-(const Vector &lhs, const Vector &rhs);
    double operator*(const Vector &lhs, const Vector &rhs);
    Vector operator*(const Vector &lhs, double rhs);
    Vector operator/(const Vector &lhs, double rhs);
    Vector operator*(double lhs, const Vector &rhs);

    bool operator==(Vector const &lhs, Vector const &rhs);
    bool operator!=(const Vector &lhs, const Vector &rhs);
}

std::ostream& operator <<(std::ostream& out, Core::Vector const& vec);


#endif /* BTree_vector_hpp */