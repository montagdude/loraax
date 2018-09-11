// Header for tripanel class

#ifndef TRIPANEL_H
#define TRIPANEL_H

#include <string>
#include <Eigen/Core>
#include "panel.h"

class Vertex;

/******************************************************************************/
//
// Tri panel class.  Derived class from panel, has 3 nodes/edges.
//
/******************************************************************************/
class TriPanel: public Panel {

    private:

    // Computing geometric quantities

    void computeCharacteristicLength ();
    void computeArea ();
    void computeNormal ();
    void computeCentroid ();
    void computeTransform ();

    public:

    // Constructor and destructor

    TriPanel ();
    virtual ~TriPanel ();

    // Add vertices and comput geometric quantities when 4 are set

    int addVertex ( Vertex * vert, bool ref_element_to_vert=true );

    // Recomputes all geometric quantities

    int recomputeGeometry ();

    // Computing source and doublet influence coefficients at a point

    double sourcePhiCoeff ( const double & x, const double & y,
                            const double & z, bool onpanel,
                            const std::string & side,
                            bool mirror_y=false ) const;
    Eigen::Vector3d sourceVCoeff ( const double & x, const double & y,
                                   const double & z, bool onpanel,
                                   const std::string & side,
                                   bool mirror_y=false ) const;
    double doubletPhiCoeff ( const double & x, const double & y,
                             const double & z, bool onpanel,
                             const std::string & side,
                             bool mirror_y=false ) const;
    Eigen::Vector3d doubletVCoeff ( const double & x, const double & y,
                                    const double & z, bool onpanel,
                                    const std::string & side,
                                    bool mirror_y=false ) const;

    // Computing induced velocity potential and induced velocity at a point

    double inducedPotential ( const double & x, const double & y,
                              const double & z, bool onpanel,
                              const std::string & side,
                              bool mirror_y=false ) const;
    Eigen::Vector3d inducedVelocity ( const double & x, const double & y,
                                      const double & z, bool onpanel,
                                      const std::string & side,
                                      bool mirror_y=false ) const;

    // Returns induced velocity due to modeling the doublet panel as a vortex
    // ring, which should give the same result as the doublet panel, but is
    // better behaved and cheaper to compute. Does not include the source
    // contribution.

    Eigen::Vector3d vortexVelocity ( const double & x, const double & y,
                                     const double & z, const double & rcore,
                                     bool mirror_y=false ) const;
};

#endif
