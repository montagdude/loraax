// Header for triface class

#ifndef TRIFACE_H
#define TRIFACE_H

#include <string>
#include "face.h"

/******************************************************************************/
//
// Tri face class.  Derived class from Face, has 3 nodes/edges.
//
/******************************************************************************/
class TriFace: public Face {

  private:

    // Computing geometric quantities

    void computeCharacteristicLength();
    void computeArea ();
    void computeNormal ();
    void computeCentroid ();
    void computeTransform ();

  public:

    // Constructor

    TriFace ();

    // Computes face geometric quantities once four nodes are set

    void setNode ( unsigned int, Node * );

    // Computing source and doublet influence coefficients at a point

    double sourcePhiCoeff ( const double &, const double &, 
                            const double &, const bool, 
                            const std::string & ) const;
    Eigen::Vector3d sourceVCoeff ( const double &, const double &,
                                   const double &, const bool,
                                   const std::string & ) const;
    double doubletPhiCoeff ( const double &, const double &, 
                             const double &, const bool,
                             const std::string & ) const;
    Eigen::Vector3d doubletVCoeff ( const double &, const double &,
                                    const double &, const bool,
                                    const std::string & ) const;

    // Computing induced velocity potential and induced velocity at a point

    double inducedPotential ( const double &, const double &,
                              const double &, const bool,
                              const std::string & ) const;
    Eigen::Vector3d inducedVelocity ( const double &, const double &,
                                      const double &, const bool,
                                      const std::string & ) const;

};

#endif
