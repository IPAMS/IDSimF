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
 PSim_simionPotentialArray.hpp

 Interface to binary electrostatic SIMION potential arrays
 (Magnetic potential arrays are not supported yet)

 ****************************/
#ifndef Particle_simulation_simion_potential_array
#define Particle_simulation_simion_potential_array


#include "Core_vector.hpp"
#include <iostream>
#include <string>
#include <array>
#include <vector>


namespace ParticleSimulation{

    using index_t = int;

    enum PASymmetry{
        CYLINDRICAL,
        PLANAR
    };

    enum PADimensionality{
        PA_2D,
        PA_3D
    };

    /**
     * Individual exception class for problems in potential arrays
     */
    class PotentialArrayException : public std::runtime_error {
    public:
        explicit PotentialArrayException (const std::string msg): std::runtime_error(msg) {}
        operator std::string() const; //Catch expects that exceptions can be cast to std::string for textual comparison
    };


    class SimionPotentialArray {

    public:

        SimionPotentialArray(std::string filename, double spatialScale = 0.001);
        SimionPotentialArray(std::string filename, double spatialScale, double potentialScale);
        SimionPotentialArray(std::string filename, Core::Vector position, double spatialScale, double potentialScale);

        //disable copies:
        SimionPotentialArray(const SimionPotentialArray&) = delete;

        double getPotential(index_t ix, index_t iy, index_t iz) const;
        double getInterpolatedPotential(double xPt, double yPt, double zPt);
        Core::Vector getField(double xPt, double yPt, double zPt);
        bool isElectrode(double xPt, double yPt, double zPt);
        bool isInside(double xPt, double yPt, double zPt);

        std::array<double, 6> getBounds();
        std::array<index_t, 3> getNumberOfGridPoints();
        std::string getHeaderString();
        void printState();

    private:

        //private methods:
        void readPa_(std::string filename);
        void readBinaryPa_(std::ifstream& inStream);
        void updateBounds_();
        std::array<double,3> transformAndCheckCoordinates_(double xPt, double yPt, double zPt) const;
        std::array<double,3> scaleAndShiftCoordinates_(double xPt, double yPt, double zPt) const;

        double interpolatedPotentialCylindrical_(double xT, double rT) const;
        double interpolatedPotentialCartesian2D_(double xT, double yT) const;
        double interpolatedPotentialCartesian3D_(double xT, double yT, double zT) const;
        double potential_(index_t ix, index_t iy, index_t iz) const;
        double rawPotential_(index_t ix, index_t iy, index_t iz) const;
        size_t linearIndex_(index_t  ix, index_t  iy, index_t  iz) const;
        bool isInside_(double x, double y, double z) const;
        bool isElectrode_(index_t ix, index_t iy, index_t iz);

        //members:
        std::string paFilename_; ///< Name of the PA file

        std::array<double,3> cornerLocation_; ///< Position of the "lower" corner of the PA in real space (to move the PA)
        double spatialScale_; ///< Geometric scale between real world coordinates and PA coordinates
        double potentialScale_; ///< Scale factor between internal potentials and real world potential

        std::vector<double> points_; ///< The raw data points of the PA in a linearized vector

        int mode_ = 0;                       ///< Mode parameter of the SIMION PA
        const double nodeDistance_ = 1.0;    ///< Distance between nodes in the PA (in internal coordinates)
        const double nodeDistanceHalf_ = nodeDistance_ / 2.0; ///< Half distance between nodes in the PA
        PASymmetry symmetry_ = CYLINDRICAL;  ///< Potential array symmetry
        PADimensionality dim_ = PA_2D;       ///< Potential array dimensionality
        double maxVoltage_ = 0.0;            ///< Maximum voltage parameter of the SIMION PA

        index_t nx_ = 0; ///< Number of nodes in x direction
        index_t ny_ = 0; ///< Number of nodes in y direction
        index_t nz_ = 0; ///< Number of nodes in z direction

        double dx_ = 0;
        double dy_ = 0;
        double dz_ = 0;

        std::array<index_t ,3> internalBoundsLower_ = {0, 0, 0};
        std::array<index_t ,3> internalBoundsUpper_ = {0, 0, 0};
        std::array<double ,6> bounds_ = {0, 0, 0, 0, 0, 0};

        int numPoints_ = 0; ///< Total number of nodes / points
        bool mirrorx_ = false; ///< Flag if mirrored in x direction
        bool mirrory_ = false; ///< Flag if mirrored in y direction
        bool mirrorz_ = false; ///< Flag if mirrored in z direction

    };
}

#endif //Particle_simulation_simion_potential_array
