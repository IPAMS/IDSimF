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

#include "Core_vector.hpp"
#include <iostream>
#include <cmath>

/**
 * Constructs vector with components
 * @param x x component
 * @param y y component
 * @param z z component
 */
Core::Vector::Vector(double x,double y, double z)
:
    x_(x),
    y_(y),
    z_(z)
{}

/**
 * Constructs vector with components passed as array
 * @param coord three element array with x,y,z components
 */
Core::Vector::Vector(const double* coord){
    x_= coord[0];
    y_= coord[1];
    z_= coord[2];
}

/**
 * Prints vector state / vector components
 */
void Core::Vector::printState() const{
    std::cout << "x="<<x_<< " y="<<y_<<" z="<<z_<<"\n";
}

// Accessors:

/**
 * x component
 */
double Core::Vector::x() const{
    return x_;
}

/**
 * y component
 */
double Core::Vector::y() const{
    return y_;
}

/**
 * z component
 */
double Core::Vector::z() const{
    return z_;
}

// Setters:

/**
 * Sets new components
 * @param x x component
 * @param y y component
 * @param z z component
 */
void Core::Vector::set(double x, double y, double z){
    x_=x;
    y_=y;
    z_=z;
}

/**
 * Set new x component
 */
void Core::Vector::x(double x){
    x_ = x;
}

/**
 * Set new y component
 */
void Core::Vector::y(double y){
    y_ = y;
}

/**
 * Set new z component
 */
void Core::Vector::z(double z){
    z_ = z;
}

// mathematical methods:

/**
 * Returns the magnitude of the vector
 */
double Core::Vector::magnitude() const{
    return sqrt(x_*x_ + y_*y_ + z_*z_);
}

/**
 * Returns the square of magnitude of the vector
 */
double Core::Vector::magnitudeSquared() const{
    return (x_*x_ + y_*y_ + z_*z_);
}

// overloaded operators:
Core::Vector Core::operator+(const Core::Vector &lhs, const Core::Vector & rhs){
    return Core::Vector(lhs.x_ + rhs.x_, lhs.y_ + rhs.y_, lhs.z_ + rhs.z_);
}

Core::Vector Core::operator-(const Core::Vector &lhs, const Core::Vector &rhs){
    return Core::Vector(lhs.x_ - rhs.x_, lhs.y_ - rhs.y_, lhs.z_ - rhs.z_);
}

double Core::operator*(const Core::Vector &lhs, const Core::Vector &rhs){
    return lhs.x_ * rhs.x_ + lhs.y_ * rhs.y_ + lhs.z_ * rhs.z_;
}

Core::Vector Core::operator*(const Core::Vector &lhs, double rhs){
    return Core::Vector(lhs.x_ * rhs, lhs.y_ * rhs, lhs.z_ * rhs);
}

Core::Vector Core::operator/(const Core::Vector &lhs, double rhs){
    return lhs * (1.0 / rhs);
}

bool Core::operator==(Core::Vector const &lhs, Core::Vector const &rhs){
    return (lhs.x_ == rhs.x_ && lhs.y_ == rhs.y_ && lhs.z_ == rhs.z_);
}

bool Core::operator!=(const Core::Vector &lhs, const Core::Vector &rhs){
    return !(lhs == rhs);
}

std::ostream& operator<< (std::ostream& os,Core::Vector const& vec)
{
    os << vec.x() << ' ' << vec.y() << ' ' << vec.z();
    return os;
}
