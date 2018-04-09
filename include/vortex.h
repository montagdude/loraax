// Header for VortexRing class

#ifndef VORTEX_H
#define VORTEX_H

#include <vector>
#include <Eigen/Core>
#include "element.h"

/******************************************************************************/
//
// Vortex class. Can be vortex ring or horseshoe vortex.
//
/******************************************************************************/
class Vortex: public Element {

  protected:

    double _gamma;

  public:

    // Constructor

    Vortex ();

    // Setting and accessing circulation

    void setCirculation ( const double & gamma );
    const double & circulation () const;

    // Induced velocity coefficient and induced velocity at a point

    virtual Eigen::Vector3d VCoeff ( const double & x, const double & y,
                                     const double & z, const double & rcore,
                                     bool mirror_y=false ) const = 0;

    virtual Eigen::Vector3d inducedVelocity ( const double & x, const double & y,
                                              const double & z,
                                              const double & rcore,
                                              bool mirror_y=false ) const = 0;
};

#endif
