// Header for panel class

#ifndef PANEL_H
#define PANEL_H

#include <vector>
#include <string>
#include <Eigen/Core>
#include "element.h"

class Vertex;

/******************************************************************************/
//
// Panel class. Computes geometric quantities, source/doublet influence
// coefficients at a point, etc.
//
/******************************************************************************/
class Panel: public Element {

  protected:

    double _sigma, _mu;            // Source and doublet strength
    double _length;                // A characteristic length
    double _area;                  // Face area
    Eigen::Vector3d _norm;         // Outward-facing normal
    Eigen::Vector3d _cen;          // Face centroid
    Eigen::Matrix3d _trans, _invtrans;
                                   // Transform from inertial frame to panel
                                   // frame and vice versa
    std::vector<double> _xtrans, _ytrans;
                                   // Panel endpoint coordinates in panel frame
    Eigen::Vector3d _vel;	   // Flow velocity at centroid

  public:

    // Constructor

    Panel ();

    // Setting and accessing source and doublet strength

    void setSourceStrength ( const double & );
    const double & sourceStrength () const;

    void setDoubletStrength ( const double & );
    const double & doubletStrength () const;

    // Returns geometric quantities

    const double & area () const;
    const Eigen::Vector3d & centroid () const;
    const Eigen::Vector3d & normal () const;

    // Setting and accessing flow velocity

    void setVelocity ( const Eigen::Vector3d & vel );
    const Eigen::Vector3d & velocity () const;

    // Virtual functions implemented in derived classes
    
    virtual double sourcePhiCoeff ( const double & x, const double & y,
                                    const double & z, bool onpanel,
                                    const std::string & side,
                                    bool mirror_y=false ) const = 0;
    virtual Eigen::Vector3d sourceVCoeff ( const double & x, const double & y,
                                           const double & z, bool onpanel,
                                           const std::string & side,
                                           bool mirror_y=false ) const = 0;
    virtual double doubletPhiCoeff ( const double & x, const double & y,
                                     const double & z, bool onpanel,
                                     const std::string & side,
                                     bool mirror_y=false ) const = 0;
    virtual Eigen::Vector3d doubletVCoeff ( const double & x, const double & y,
                                           const double & z, bool onpanel,
                                           const std::string & side,
                                           bool mirror_y=false ) const = 0;
    virtual double inducedPotential ( const double & x, const double & y,
                                      const double & z, bool onpanel,
                                      const std::string & side,
                                      bool mirror_y=false ) const = 0;
    virtual Eigen::Vector3d inducedVelocity ( const double & x,
                                           const double & y, const double & z,
                                           bool onpanel,
                                           const std::string & side,
                                           bool mirror_y=false ) const = 0;
};

#endif
