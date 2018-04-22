// Header for panel class

#ifndef PANEL_H
#define PANEL_H

#include <vector>
#include <string>
#include <Eigen/Dense>
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

    Panel * _right, * _left, * _front, * _back;
				   // Neighbor panels
    Eigen::Matrix3d _jac;	   // Grid metrics jacobian matrix
    Eigen::PartialPivLU<Eigen::Matrix3d> _lu;
				   // Factorization of jacobian matrix 

    // Computing geometric quantities

    virtual void computeCharacteristicLength () = 0;
    virtual void computeArea () = 0;
    virtual void computeNormal () = 0;
    virtual void computeCentroid () = 0;
    virtual void computeTransform () = 0;

  public:

    // Constructor

    Panel ();

    // Recomputes all geometric quantities

    virtual int recomputeGeometry () = 0;

    // Set/access neighbor panels

    void setRightNeighbor ( Panel * right ); 
    void setLeftNeighbor ( Panel * left ); 
    void setFrontNeighbor ( Panel * front ); 
    void setBackNeighbor ( Panel * back ); 

    Panel * rightNeighbor ();
    Panel * leftNeighbor ();
    Panel * frontNeighbor ();
    Panel * backNeighbor ();

    // Compute grid transformation (must be done after setting neighbors)

    int computeGridTransformation ();

    // Setting and accessing source and doublet strength

    void setSourceStrength ( const double & );
    const double & sourceStrength () const;

    void setDoubletStrength ( const double & );
    const double & doubletStrength () const;

    // Returns geometric quantities

    const double & area () const;
    const Eigen::Vector3d & centroid () const;
    const Eigen::Vector3d & normal () const;

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

    // Compute or access surface velocity. Must set neighbors and compute grid
    // transformation before computing velocity.

    void computeVelocity ( const Eigen::Vector3d & uinfvec );
    const Eigen::Vector3d & velocity () const;
};

#endif
