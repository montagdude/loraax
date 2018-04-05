// Header for quadpanel class

#ifndef QUADPANEL_H
#define QUADPANEL_H

#include <string>
#include <Eigen/Core>
#include "panel.h"

class Vertex;

/******************************************************************************/
//
// Quad panel class.  Derived class from Panel, has 4 endpoints.
//
/******************************************************************************/
class QuadPanel: public Panel {

  private:

    // Computing geometric quantities

    void computeCharacteristicLength();
    void computeArea ();
    void computeNormal ();
    void computeCentroid ();
    void computeTransform ();

  public:

    // Constructor

    QuadPanel ();

    // Add vertices and compute geometric quantities when 4 are set

    int addVertex ( Vertex * vert );

    // Computing source and doublet influence coefficients at a point

    double sourcePhiCoeff ( const double & x, const double & y,
                            const double & z, const bool onpanel,
                            const std::string & side ) const;
    Eigen::Vector3d sourceVCoeff ( const double & x, const double & y,
                                   const double & z, const bool onpanel,
                                   const std::string & side ) const;
    double doubletPhiCoeff ( const double & x, const double & y,
                             const double & z, const bool onpanel,
                             const std::string & side ) const;
    Eigen::Vector3d doubletVCoeff ( const double & x, const double & y,
                                    const double & z, const bool onpanel,
                                    const std::string & side ) const;

    // Computing induced velocity potential and induced velocity at a point

    double inducedPotential ( const double & x, const double & y,
                              const double & z, const bool onpanel,
                              const std::string & side ) const;
    Eigen::Vector3d inducedVelocity ( const double & x, const double & y,
                                      const double & z, const bool onpanel,
                                      const std::string & side ) const;

};

#endif
