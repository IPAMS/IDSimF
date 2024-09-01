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

#include "Core_debug.hpp"
#include "Core_utils.hpp"
#include "PSim_simionPotentialArray.hpp"
#include "spdlog/spdlog.h"
#include <sstream>
#include <cmath>

ParticleSimulation::PotentialArrayException::operator std::string() const {
    return this->what();
}

/**
 * Default constructor: Reads a SIMION PA file. The geometric scale is 0.001
 * since the usual length unit in SIMION PA files is mm.
 *
 * @param filename The filename of the SIMION PA file to read
 * @param spatialScale A geometric scale factor to scale the potential array
 * (default = 0.001 which is the conversion factor from mm, the length unit of SIMION, to m)
 */
ParticleSimulation::SimionPotentialArray::SimionPotentialArray(std::string filename, double spatialScale)
        :cornerLocation_({0.0, 0.0, 0.0}), spatialScale_(spatialScale), potentialScale_(1.0) {
    readPa_(filename);
}

ParticleSimulation::SimionPotentialArray::SimionPotentialArray(std::string filename,
                                                               double spatialScale, double potentialScale)
        :cornerLocation_({0.0, 0.0, 0.0}), spatialScale_(spatialScale), potentialScale_(potentialScale) {
    readPa_(filename);
}


ParticleSimulation::SimionPotentialArray::SimionPotentialArray(std::string filename, Core::Vector position,
                                                               double spatialScale, double potentialScale)
        :cornerLocation_({position.x(),position.y(),position.z()}), spatialScale_(spatialScale), potentialScale_(potentialScale) {
    readPa_(filename);
}

void ParticleSimulation::SimionPotentialArray::setCornerPosition(Core::Vector newCornerPosition) {
    cornerLocation_ = {newCornerPosition.x(),newCornerPosition.y(),newCornerPosition.z()};
    updateBounds_();
}


double ParticleSimulation::SimionPotentialArray::getPotential(index_t ix, index_t iy, index_t iz) const {

    if (Core::safetyGuards) {
        if (
                ix<0 || ix>=nx_ ||
                iy<0 || iy>=ny_ ||
                iz<0 || iz>=nz_) {
            std::stringstream ss;
            ss << "Index " << ix << " " << iy << " " << iz << " is not in the potential array";
            throw (ParticleSimulation::PotentialArrayException(ss.str()));
        }
    }

    double pot = potential_(ix, iy, iz);
    // electrodes are defined by potentials larger than maxVoltage_ in simion PAs
    if (pot > maxVoltage_){
        pot -= 2*maxVoltage_;
    }
    return pot;
}

double ParticleSimulation::SimionPotentialArray::getInterpolatedPotential(double xPt, double yPt, double zPt) const{

    double pot = 0.0;
    if (symmetry_ == CYLINDRICAL){

        //calculate normalized grid coordinates:
        std::array<double,3> tc= transformAndCheckCoordinates_(xPt, yPt, zPt);
        double xT = tc[0];
        double rT = tc[1];

        pot = interpolatedPotentialCylindrical_(xT, rT);

    } else  { //symmetry: Planar

        //calculate normalized grid coordinates:
        std::array<double,3> tc= transformAndCheckCoordinates_(xPt, yPt, zPt);
        double xT = tc[0];
        double yT = tc[1];
        double zT = tc[2];

        if (dim_ == PA_2D){
            pot = interpolatedPotentialCartesian2D_(xT, yT);
        }
        else {
            pot = interpolatedPotentialCartesian3D_(xT, yT, zT);
        }
    }
    return pot;
}


Core::Vector ParticleSimulation::SimionPotentialArray::getField(double xPt, double yPt, double zPt) const{

     Core::Vector field;
    if (symmetry_ == CYLINDRICAL){
        //calculate normalized grid coordinates:
        std::array<double,3> tc= scaleAndShiftCoordinates_(xPt, yPt, zPt);
        double xT = tc[0];
        double yT = tc[1];
        double zT = tc[2];
        double rT = sqrt(yT*yT + zT*zT);

        if (Core::safetyGuards && (xT < internalBoundsLower_[0] ||
                                   xT > internalBoundsUpper_[0] ||
                                   rT > internalBoundsUpper_[1] ) ){
            std::stringstream ss;
            ss << "Point " << xPt << " " << yPt << " " << zPt << " is not in the cylindrical potential array";
            throw (ParticleSimulation::PotentialArrayException(ss.str()));
        }

        // Handle grid boundaries (use value on boundary if difference quotient reaches over the boundary)
        double xL = xT - nodeDistanceHalf_;
        if (xL < internalBoundsLower_[0]) xL = internalBoundsLower_[0];

        double xU = xT + nodeDistanceHalf_;
        if (xU > internalBoundsUpper_[0]) xU = internalBoundsUpper_[0];

        double rL = rT - nodeDistanceHalf_;
        if (rL < internalBoundsLower_[1]) rL = internalBoundsLower_[1];

        double rU = rT + nodeDistanceHalf_;
        if (rU > internalBoundsUpper_[1]) rU = internalBoundsUpper_[1];

        double pL0 = interpolatedPotentialCylindrical_(xL, rT);
        double pU0 = interpolatedPotentialCylindrical_(xU, rT);
        double p0L = interpolatedPotentialCylindrical_(xT, rL);
        double p0U = interpolatedPotentialCylindrical_(xT, rU);

        double Ex = (pL0 - pU0) / (xU-xL);
        double Er = (p0L - p0U) / (rU-rL);

        // the radial field points outwards radially, since we can scale its cartesian components
        // with the normalized cartesian components of the point

        if (Core::isDoubleUnequal(rT, 0.0)){
            double Ey = Er * yT / rT;
            double Ez = Er * zT / rT;
            field = {Ex, Ey, Ez};
        }
        else {
            // at exactly the center line: there is no radial field
            field = {Ex, 0.0, 0.0};
        }
    }
    else {
        std::array<double,3> tc= transformAndCheckCoordinates_(xPt, yPt, zPt);
        double xT = tc[0];
        double yT = tc[1];
        double zT = tc[2];

        // Handle grid boundaries (use value on boundary if difference quotient reaches over the boundary)
        double xL = xT - nodeDistanceHalf_;
        if (xL < internalBoundsLower_[0]) xL = internalBoundsLower_[0];

        double xU = xT + nodeDistanceHalf_;
        if (xU > internalBoundsUpper_[0]) xU = internalBoundsUpper_[0];

        double yL = yT - nodeDistanceHalf_;
        if (yL < internalBoundsLower_[1]) yL = internalBoundsLower_[1];

        double yU = yT + nodeDistanceHalf_;
        if (yU > internalBoundsUpper_[1]) yU = internalBoundsUpper_[1];

        if (dim_ == PA_2D){

            double pL0 = interpolatedPotentialCartesian2D_(xL, yT);
            double pU0 = interpolatedPotentialCartesian2D_(xU, yT);
            double p0L = interpolatedPotentialCartesian2D_(xT, yL);
            double p0U = interpolatedPotentialCartesian2D_(xT, yU);

            double Ex = (pL0 - pU0) / (xU-xL);
            double Ey = (p0L - p0U) / (yU-yL);

            field = {Ex, Ey, 0.0};
        }
        else {
            double zL = zT - nodeDistanceHalf_;
            if (zL < internalBoundsLower_[2]) zL = internalBoundsLower_[2];

            double zU = zT + nodeDistanceHalf_;
            if (zU > internalBoundsUpper_[2]) zU = internalBoundsUpper_[2];

            double pL00 = interpolatedPotentialCartesian3D_(xL, yT, zT);
            double pU00 = interpolatedPotentialCartesian3D_(xU, yT, zT);
            double p0L0 = interpolatedPotentialCartesian3D_(xT, yL, zT);
            double p0U0 = interpolatedPotentialCartesian3D_(xT, yU, zT);
            double p00L = interpolatedPotentialCartesian3D_(xT, yT, zL);
            double p00U = interpolatedPotentialCartesian3D_(xT, yT, zU);

            double Ex = (pL00 - pU00) / (xU-xL);
            double Ey = (p0L0 - p0U0) / (yU-yL);
            double Ez = (p00L - p00U) / (zU-zL);

            field = {Ex, Ey, Ez};
        }
    }

    // The field is expected to be returned in V/m and has therefore to be
    // geometrically scaled
    return(field / spatialScale_);
}


bool ParticleSimulation::SimionPotentialArray::isElectrode(double xPt, double yPt, double zPt) const{
    std::array<double,3> posTransformed = transformAndCheckCoordinates_(xPt, yPt, zPt);
    double xT = posTransformed[0];
    double yT = posTransformed[1];
    double zT = posTransformed[2];

    if (mirrorx_ &&  xT < 0.0)
        xT = -xT;
    if (mirrory_ && yT < 0.0)
        yT = -yT;
    if (mirrorz_ && zT < 0.0)
        zT = -zT;

    index_t xNode = static_cast<index_t>(xT);
    index_t yNode = static_cast<index_t>(yT);
    index_t zNode = static_cast<index_t>(zT);


    // as soon as one coordinate is not exactly on a node coordinate,
    // we have to check if the surrounding nodes are part of an electrode,
    // which makes the current point part of an electrode
    // (xT and xNode for example would be exactly equal if the transformed coordinate is exactly on a node)

    if (xT > xNode || yT > yNode || zT > zNode){

        // as soon one point of the surrounding points is NOT an electrode point:
        // The current point is not in an electrode
        if (dim_ == PA_2D){

            // since "point is in an electrode" is defined by the surrounding grid nodes being electrode
            // and we look to the right / upwards when being exactly on an edge between nodes, the right / upper
            // boundary of the grid has to be considered separately and is defined as "not electrode"
            //
            if (xNode == nx_-1 || yNode == ny_-1){
                return false;
            }
            else{
                return !(!isElectrode_(xNode,   yNode,   zNode) ||
                         !isElectrode_(xNode+1, yNode,   zNode) ||
                         !isElectrode_(xNode,   yNode+1, zNode) ||
                         !isElectrode_(xNode+1, yNode+1, zNode));
            }
        }
        else {

            // Consider upper boundaries:
            if (xNode == nx_-1 || yNode == ny_-1 || zNode == nz_-1){
                return false;
            }
            else {
                return !(!isElectrode_(xNode,  yNode,   zNode) ||
                         !isElectrode_(xNode+1, yNode,   zNode) ||
                         !isElectrode_(xNode,   yNode+1, zNode) ||
                         !isElectrode_(xNode+1, yNode+1, zNode) ||
                         !isElectrode_(xNode,   yNode,   zNode+1) ||
                         !isElectrode_(xNode+1, yNode,   zNode+1) ||
                         !isElectrode_(xNode,   yNode+1, zNode+1) ||
                         !isElectrode_(xNode+1, yNode+1, zNode+1)
                );
            }
        }
    }
    else {
        return isElectrode_(xNode, yNode, zNode);
    }
}


bool ParticleSimulation::SimionPotentialArray::isInside(double xPt, double yPt, double zPt) const{
    std::array<double,3> posTransformed = scaleAndShiftCoordinates_(xPt, yPt, zPt);
    return isInside_(posTransformed[0], posTransformed[1], posTransformed[2]);
}

std::array<double,6> ParticleSimulation::SimionPotentialArray::getBounds() const{
    return bounds_;
}

std::array<ParticleSimulation::index_t,3> ParticleSimulation::SimionPotentialArray::getNumberOfGridPoints() const{
    return {nx_,ny_,nz_};
}

std::string  ParticleSimulation::SimionPotentialArray::getHeaderString() const{
    return "not yet implemented";
}

void ParticleSimulation::SimionPotentialArray::readPa_(std::string filename) {
    auto logger = spdlog::get("logger");
    if (logger){
        logger->info("Reading PA: {}", filename);
    }
    paFilename_ = filename;
    std::ifstream inStream(filename, std::ios::binary);
    readBinaryPa_(inStream);
    updateBounds_();
}

void ParticleSimulation::SimionPotentialArray::readBinaryPa_(std::ifstream &inStream) {
    //
    if (!inStream.read(reinterpret_cast<char*> (&mode_), sizeof(mode_))){
        throw (ParticleSimulation::PotentialArrayException("File empty"));
    };

    int symmetryBuf = 0;
    inStream.read( reinterpret_cast<char*>(&symmetryBuf), sizeof(symmetryBuf));
    inStream.read( reinterpret_cast<char*>(&maxVoltage_), sizeof(maxVoltage_));

    unsigned rawMirror = 0;

    inStream.read( reinterpret_cast<char*>(&nx_), sizeof(nx_));
    inStream.read( reinterpret_cast<char*>(&ny_), sizeof(ny_));
    inStream.read( reinterpret_cast<char*>(&nz_), sizeof(nz_));
    inStream.read( reinterpret_cast<char*>(&rawMirror), sizeof(unsigned));

    //FIXME: use bitmasks here...
    mirrorx_ = (rawMirror & 1) == 1;
    mirrory_ = ((rawMirror >> 1) & 1) == 1;
    mirrorz_ = ((rawMirror >> 2) & 1) == 1;

    bool electrostatic = !(((rawMirror >> 3) & 1) == 1);
    //----------
    if (!electrostatic){
        throw (ParticleSimulation::PotentialArrayException(
                "Simion PA file is not electrostatic. Currently only electrostatic PAs are supported"));
    }

    //unsigned ng = (rawMirror >> 4) & ((1 << 17)-1);

    // set symmetry and dimensionality:
    if (symmetryBuf == 0){
        symmetry_ = CYLINDRICAL;
        dim_ = PA_2D;

        if (mirrorx_){
            internalBoundsLower_ = {-(nx_ -1), 0, 0};
        } else{
            internalBoundsLower_ = {0, 0, 0};
        }
        internalBoundsUpper_ = {nx_-1, ny_-1, 0};

    } else {
        symmetry_ = PLANAR;

        index_t lbx = 0;
        index_t lby = 0;
        index_t lbz = 0;

        if (mirrorx_){
            lbx = -(nx_ -1);
        }
        if (mirrory_){
            lby = -(ny_ -1);
        }
        if (mirrorz_){
            lbz = -(nz_ -1);
        }

        internalBoundsLower_ = {lbx, lby, lbz};
        internalBoundsUpper_ = {nx_-1, ny_-1, nz_-1};

        if (nz_ == 1){
            dim_ = PA_2D;
        } else {
            dim_ = PA_3D;
        }
    }

    // read internal scaling ...
    if (mode_ <= -2){
        inStream.read( reinterpret_cast<char*>(&dx_), sizeof(dx_));
        inStream.read( reinterpret_cast<char*>(&dy_), sizeof(dy_));
        inStream.read( reinterpret_cast<char*>(&dz_), sizeof(dz_));
    } else {
        dx_ = 1.0;
        dy_ = 1.0;
        dz_ = 1.0;
    }

    numPoints_ = static_cast<unsigned int>(nx_ * ny_ * nz_);
    points_.resize(numPoints_);
    inStream.read( reinterpret_cast<char*>(points_.data()), static_cast<std::streamsize >(sizeof(double)*numPoints_));
}


/**
 * Updates bounds, which are the extent of the PA in the spatial directions.
 */
void ParticleSimulation::SimionPotentialArray::updateBounds_() {
    if (symmetry_ == CYLINDRICAL){
        double xMin = cornerLocation_[0] + internalBoundsLower_[0] * spatialScale_;
        double xMax = cornerLocation_[0] + internalBoundsUpper_[0] * spatialScale_;
        double yMin = cornerLocation_[1] - internalBoundsUpper_[1] * spatialScale_;
        double yMax = cornerLocation_[1] + internalBoundsUpper_[1] * spatialScale_;
        double zMin = cornerLocation_[2] - internalBoundsUpper_[1] * spatialScale_;
        double zMax = cornerLocation_[2] + internalBoundsUpper_[1] * spatialScale_;

        bounds_ = {xMin, xMax, yMin, yMax, zMin, zMax};
    }
    else{
        double xMin = cornerLocation_[0] + internalBoundsLower_[0] * spatialScale_;
        double xMax = cornerLocation_[0] + internalBoundsUpper_[0] * spatialScale_;
        double yMin = cornerLocation_[1] + internalBoundsLower_[1] * spatialScale_;
        double yMax = cornerLocation_[1] + internalBoundsUpper_[1] * spatialScale_;
        double zMin = cornerLocation_[2] + internalBoundsLower_[2] * spatialScale_;
        double zMax = cornerLocation_[2] + internalBoundsUpper_[2] * spatialScale_;

        bounds_ = {xMin, xMax, yMin, yMax, zMin, zMax};
    }

}


/**
 * Calculates the transformed / mapped grid coordinates (including translation, scaling and mirroring) from a
 * world coordinate
 */
std::array<double,3> ParticleSimulation::SimionPotentialArray::transformAndCheckCoordinates_(double xPt, double yPt, double zPt) const{

    double xT = (xPt-cornerLocation_[0]) / spatialScale_;
    double yT = (yPt-cornerLocation_[1]) / spatialScale_;
    double zT = (zPt-cornerLocation_[2]) / spatialScale_;

    if (symmetry_ == CYLINDRICAL){

        double r = sqrt(yT*yT + zT*zT);

        if (Core::safetyGuards && (xT < internalBoundsLower_[0] ||
                                   xT > internalBoundsUpper_[0] ||
                                   r  > internalBoundsUpper_[1] ) ){
            std::stringstream ss;
            ss << "Point " << xPt << " " << yPt << " " << zPt << " is not in the cylindrical potential array";
            throw (ParticleSimulation::PotentialArrayException(ss.str()));
        }

        if (mirrorx_ &&  xT < 0.0)
            xT = -xT;
        return {xT, r, 0.0};
    }
    else{

        if (Core::safetyGuards &&
                (xT < internalBoundsLower_[0] || xT > internalBoundsUpper_[0] ||
                 yT < internalBoundsLower_[1] || yT > internalBoundsUpper_[1] ||
                 zT < internalBoundsLower_[2] || zT > internalBoundsUpper_[2] )) {
            std::stringstream ss;
            ss << "Point " << xPt << " " << yPt << " " << zPt << " is not in the planar potential array";
            throw (ParticleSimulation::PotentialArrayException(ss.str()));
        }

        return {xT,yT,zT};
    }
}

/**
 * Calculate translated and scaled coordinates from world coordinate
 */
std::array<double,3> ParticleSimulation::SimionPotentialArray::scaleAndShiftCoordinates_(double xPt, double yPt, double zPt) const{

    double xt = (xPt-cornerLocation_[0]) / spatialScale_;
    double yt = (yPt-cornerLocation_[1]) / spatialScale_;
    double zt = (zPt-cornerLocation_[2]) / spatialScale_;

    return  {xt,yt,zt};
}

/**
 * Calculates the interpolated potential in cylindrical coordinates
 */
double ParticleSimulation::SimionPotentialArray::interpolatedPotentialCylindrical_(double xT, double rT) const {
    // in cylindrical coordinates
    // mirror in positive internal coordinates if necessary:
    if (mirrorx_ &&  xT < 0.0){
        xT = -xT;
    }

    index_t xNode = static_cast<index_t>(xT);
    index_t rNode = static_cast<index_t>(rT);

    double xWeight = xT - xNode;
    double rWeight = rT - rNode;

    // get the values on the 4 nodes surrounding the point we want to get the potential for:
    // if the weight in one direction is 0.0, we are exactly on a node position ni,
    // thus the contribution of the ni+1th node is zero anyways and we don't need the
    // potential from that node. Note that this also protects from accessing non existing out of bounds nodes
    // if we are exactly on the last node position in a direction

    //note that "potential_" makes NO bounds check (getPotential does)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
    double p00 = potential_(xNode, rNode, 0);
    double p10 = xWeight != 0.0 ? potential_(xNode+1, rNode, 0) : 0.0;
    double p01 = rWeight != 0.0 ? potential_(xNode, rNode+1, 0) : 0.0;
    double p11 = (xWeight != 0.0 && rWeight != 0.0) ? potential_(xNode+1, rNode+1, 0) : 0.0;
#pragma GCC diagnostic pop

    double pot =
            (1-xWeight) * (1-rWeight) * p00 +
               xWeight  * (1-rWeight) * p10 +
            (1-xWeight) * rWeight     * p01 +
               xWeight  * rWeight     * p11;

    return pot;
}

/**
 * Calculates an interpolated potential in a two dimensional cartesian PA
 * from the internal x,y,z coordinates xT, yT, zT
 */
double ParticleSimulation::SimionPotentialArray::interpolatedPotentialCartesian2D_(double xT, double yT) const {

    if (mirrorx_ &&  xT < 0.0)
        xT = -xT;
    if (mirrory_ && yT < 0.0)
        yT = -yT;

    index_t xNode = static_cast<index_t>(xT);
    index_t yNode = static_cast<index_t>(yT);

    double xWeight = xT - xNode;
    double yWeight = yT - yNode;

    // two dimensional interpolation is technically exactly the same as in cylindrical coordinates,
    // see the comment there
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
    double p00 = potential_(xNode, yNode, 0);
    double p10 = xWeight != 0.0 ? potential_(xNode+1, yNode, 0) : 0.0;
    double p01 = yWeight != 0.0 ? potential_(xNode, yNode+1, 0) : 0.0;
    double p11 = (xWeight != 0.0 && yWeight != 0.0) ? potential_(xNode+1, yNode+1, 0) : 0.0;
#pragma GCC diagnostic pop

    double pot =
                 (1-xWeight) * (1-yWeight) * p00 +
                  xWeight    * (1-yWeight) * p10 +
                 (1-xWeight) *  yWeight    * p01 +
                  xWeight    * yWeight     * p11;

    return pot;
}

/**
 * Calculates an interpolated potential in a three dimensional cartesian PA
 * from the internal x,y,z coordinates xT, yT, zT
 */
double ParticleSimulation::SimionPotentialArray::interpolatedPotentialCartesian3D_(double xT, double yT, double zT) const {

    if (mirrorx_ &&  xT < 0.0)
        xT = -xT;
    if (mirrory_ && yT < 0.0)
        yT = -yT;
    if (mirrorz_ && zT < 0.0)
        zT = -zT;

    index_t xNode = static_cast<index_t>(xT);
    index_t yNode = static_cast<index_t>(yT);
    index_t zNode = static_cast<index_t>(zT);

    double xWeight = xT - xNode;
    double yWeight = yT - yNode;
    double zWeight = zT - zNode;

    // get the values on the 4 nodes with lower z index surrounding the point we want to get the potential for
    // (see comment above)

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
    double p000 = potential_(xNode, yNode, zNode);
    double p100 = xWeight != 0.0 ? potential_(xNode+1, yNode, zNode) : 0.0;
    double p010 = yWeight != 0.0 ? potential_(xNode, yNode+1, zNode) : 0.0;
    double p110 = (xWeight != 0.0 && yWeight != 0.0) ? potential_(xNode+1, yNode+1, zNode) : 0.0;
#pragma GCC diagnostic pop

    double pot = (1-xWeight) * (1-yWeight) * (1-zWeight) * p000 +
                 xWeight     * (1-yWeight) * (1-zWeight) * p100 +
                 (1-xWeight) * yWeight     * (1-zWeight) * p010 +
                 xWeight     * yWeight     * (1-zWeight) * p110;

    // the part with the upper z index (zNode+1) only has to be calculated if we are not exactly on a
    // node position in z direction. If we are exactly on a node in z direction, the contribution for
    // all nodes with zNode+1 is zero, since the weighting in z direction would be zWeight which is 0.0
    // in this case.
    if (zWeight > 0.0){
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
        double p001 = potential_(xNode, yNode, zNode+1);
        double p101 = xWeight != 0.0 ? potential_(xNode+1, yNode, zNode+1) : 0.0;
        double p011 = yWeight != 0.0 ? potential_(xNode, yNode+1, zNode+1) : 0.0;
        double p111 = (xWeight != 0.0 && yWeight != 0.0) ? potential_(xNode+1, yNode+1, zNode+1) : 0.0;
#pragma GCC diagnostic pop

        pot +=  (1-xWeight) * (1-yWeight) * zWeight * p001 +
                xWeight     * (1-yWeight) * zWeight * p101 +
                (1-xWeight) * yWeight     * zWeight * p011 +
                xWeight     * yWeight     * zWeight * p111;
    }

    return pot;
}

double ParticleSimulation::SimionPotentialArray::potential_(index_t ix, index_t iy, index_t iz) const {
    double pot = rawPotential_(ix, iy, iz);
    // electrodes are defined by potentials larger than maxVoltage_ in simion PAs
    if (pot > maxVoltage_){
        pot -= 2*maxVoltage_;
    }
    return pot * potentialScale_;
}

double ParticleSimulation::SimionPotentialArray::rawPotential_(index_t ix, index_t iy, index_t iz) const{
    return points_[static_cast<std::size_t>(linearIndex_(ix,iy,iz))];
}

ParticleSimulation::index_t ParticleSimulation::SimionPotentialArray::linearIndex_(index_t ix, index_t iy, index_t iz) const{
    // This guard has a strong runtime penalty, therefore, the methods reachable from user side
    // has to check that no boundary violation happens...
    /*if (Core::safetyGuards) {
        if (
                ix<0 || ix>=nx_ ||
                iy<0 || iy>=ny_ ||
                iz<0 || iz>=nz_) {
            std::stringstream ss;
            ss << "Index " << ix << " " << iy << " " << iz << " is not in the potential array";
            throw (ParticleSimulation::PotentialArrayException(ss.str()));
        }
    }*/
    return (iz * ny_ + iy) * nx_ + ix;
}

/**
 * Internal isInside method (without coordinate transformation)
 */
bool ParticleSimulation::SimionPotentialArray::isInside_(double x, double y, double z) const{
    bool xIn = false;

    if (x < 0.0 && mirrorx_)
        x = -x;
    xIn = (x >= 0.0 && x<= nx_ - 1);

    if (symmetry_ == CYLINDRICAL){

        double r = sqrt(y*y + z*z);
        bool rIn = (r <= ny_ - 1);

        return (xIn && rIn);
    }
    else { // planar

        bool yIn = false;
        bool zIn = false;

        if (y < 0.0 && mirrory_)
            y = -y;
        yIn = (y >= 0.0 && y<= ny_ - 1);

        if (z < 0.0 && mirrorz_)
            z = -z;
        zIn = (z >= 0.0 && z<= nz_ - 1);

        return (xIn && yIn && zIn);
    }
}

bool ParticleSimulation::SimionPotentialArray::isElectrode_(index_t ix, index_t iy, index_t iz) const{
    return rawPotential_(ix, iy, iz) > maxVoltage_;
}

void ParticleSimulation::SimionPotentialArray::printState() const{

    std::cout << "SIMION PA ---------------------------------"<<std::endl;
    std::cout << "corner position: " << cornerLocation_[0] << " " << cornerLocation_[1] << " " << cornerLocation_[2] << " scale:" << spatialScale_ << std::endl;
    std::string symmetryString = (symmetry_ == CYLINDRICAL) ? "cylindrical " : "planar";

    std::cout << "mode: "<<mode_ << " symmetry: "<< symmetry_ <<" ("<<symmetryString<<")"<< " maxVoltage: "<<maxVoltage_ <<std::endl;
    std::cout << "nx: "<<nx_ << " ny: "<<ny_ << " nz: "<<nz_ <<" npoints: "<< numPoints_<<std::endl;
    std::cout << "mx: "<<mirrorx_ << " my: "<<mirrory_<< " mz: "<<mirrorz_ <<std::endl;
    std::cout << "bounds: ";
    for(auto &bound: getBounds()){
        std::cout << " "<<bound;
    }
    std::cout << std::endl;
    //std::cout << "electrostatic: "<< electrostatic << " ng:"<< ng << std::endl;

    for (std::size_t i=0; i< 10; ++i){
        std::cout << " "<< points_[i];
    }
    std::cout << std::endl;
}
