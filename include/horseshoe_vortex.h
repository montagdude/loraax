// Header for HorseshoeVortex class

#ifndef HORSESHOEVORTEX_H
#define HORSESHOEVORTEX_H

#include <vector>
#include <Eigen/Core>
#include "vortex.h"

class Vertex;

/******************************************************************************/
//
// Horseshoe vortex class
//
/******************************************************************************/
class HorseshoeVortex: public Vortex {

  public:

    // Constructor

    HorseshoeVortex ();

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
