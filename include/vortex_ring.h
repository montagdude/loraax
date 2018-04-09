// Header for VortexRing class

#ifndef VORTEXRING_H
#define VORTEXRING_H

#include <vector>
#include <Eigen/Core>
#include "vortex.h"

class Vertex;

/******************************************************************************/
//
// Vortex ring class
//
/******************************************************************************/
class VortexRing: public Vortex {

  public:

    // Constructor

    VortexRing ();

    // Set vertices

    int addVertex ( Vertex * vert );

    // Induced velocity coefficient and induced velocity at a point

    Eigen::Vector3d VCoeff ( const double & x, const double & y,
                             const double & z, const double & rcore,
                             bool mirror_y=false ) const;

    Eigen::Vector3d inducedVelocity ( const double & x, const double & y,
                                      const double & z,
                                      const double & rcore,
                                      bool mirror_y=false ) const;
};

#endif
