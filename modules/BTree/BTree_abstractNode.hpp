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
 BTree_node.hpp

 Abstract base class for Nodes in Barnes-Hut Tree
 ****************************/

#ifndef BTree_abstract_node_hpp
#define BTree_abstract_node_hpp


#include "Core_vector.hpp"
#include "Core_particle.hpp"
#include <ostream>

namespace BTree {

    class Particle;

    /**
     * Abstract base class for Barnes-Hut Tree nodes
     */
    class AbstractNode {

    public:

        /**
         * Sub-Octants of a node (e.g. SWT = south west top, NEB = north east bottom)
         */
        enum Octant {SWT,NWT,SET,NET,SWB,NWB,SEB,NEB};

        // Constructor / Destructor:
        AbstractNode(Core::Vector min, Core::Vector max);
        virtual ~AbstractNode() = default;

        // Copy Constructors:
        AbstractNode(const AbstractNode& that) = default;
        //AbstractNode& operator=(const AbstractNode& that);

        // Static methods:
        static int getNumberOfNodes();
        static Core::Vector calculateElectricField(const Core::Vector &r1, const Core::Vector &r2, double charge2);

        // Accessors:
        [[nodiscard]] std::size_t getNumberOfParticles() const;
        [[nodiscard]] Core::Particle* getParticle() const;
        [[nodiscard]] Core::Vector getCenter() const;
        [[nodiscard]] Core::Vector getMax() const;
        [[nodiscard]] Core::Vector getMin() const;
        [[nodiscard]] double getCharge() const;
        [[nodiscard]] Core::Vector getCenterOfCharge() const;
        [[nodiscard]] Octant getOctant(Core::Vector location) const;

        // Virtual member methods:
        virtual void initAsRoot() = 0;
        [[nodiscard]] virtual bool isRoot() const = 0;
        virtual void updateParents() = 0;
        virtual void updateSelf() = 0;
        virtual void removeMyselfFromTree() = 0;
        virtual void insertParticle(Core::Particle* particle) = 0;
        virtual void computeChargeDistributionRecursive() = 0;

        //export methods:
        [[nodiscard]] virtual std::string toString() const = 0;

    protected:
        static inline double theta = 0.9; ///< Theta (distance criterion) factor
        static inline int nNodes_ = 0;  ///< Number of total nodes existing

        //member fields:
        double charge_ = 0.0; ///< Total charge in the node (and the subnodes)
        Core::Vector centerOfCharge_ = Core::Vector(0,0,0); ///< The center of charge in the node
        std::size_t numP_ = 0; ///< Total number of particles

        Core::Vector min_ = Core::Vector(0,0,0); ///< Minimal corner (xlo,ylo,zlo corner) of the node
        Core::Vector max_ = Core::Vector(0,0,0); ///< Maximal corner (xhi,yhi,zhi corner) of the node
        Core::Vector center_ = Core::Vector(0,0,0); ///< Geometric center of the node
        Core::Particle* particle_ = nullptr; ///< A particle potentially located in the node
    };
}

#endif //BTree_abstract_node_hpp
