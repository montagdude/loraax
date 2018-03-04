// Header for quadface class

#ifndef QUADFACE_H
#define QUADFACE_H

#include <string>
#include "face.h"

/******************************************************************************/
//
// Quad face class.  Derived class from Face, has 4 endpoints.
//
/******************************************************************************/
class QuadFace: public Face {

  private:

    // Computing geometric quantities

    void computeCharacteristicLength();
    void computeArea ();
    void computeNormal ();
    void computeCentroid ();
    void computeTransform ();

  public:

    // Constructor

    QuadFace ();

    // Sets endpoints and computes face geometric quantities

    void setEndpoints( const std::vector<double> & x,
                       const std::vector<double> & y,
                       const std::vector<double> & z );

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
