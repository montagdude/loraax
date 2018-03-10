// Header for face class

#ifndef FACE_H
#define FACE_H

#include <vector>
#include <string>
#include <Eigen/Core>

class Vertex;

/******************************************************************************/
//
// Face class. Computes geometric quantities, source/doublet influence
// coefficients at a point, etc.
//
/******************************************************************************/
class Face {

  protected:

    int _idx;		   	   // Face identifier index
    double _sigma, _mu;            // Source and doublet strength
    double _length;                // A characteristic length
    double _area;                  // Face area
    Eigen::Vector3d _norm;         // Outward-facing normal
    Eigen::Vector3d _cen;          // Face centroid
    Eigen::Matrix3d _trans, _invtrans;
                                   // Transform from inertial frame to panel
                                   // frame and vice versa
    unsigned int _currverts;
    std::vector<Vertex *> _verts;  // Panel endpoint vertices
    std::vector<double> _xtrans, _ytrans;
                                   // Panel endpoint coordinates in panel frame

  public:

    // Constructor

    Face ();

    // Set / access index

    void setIdx ( int idx );
    int idx () const;

    // Set / access vertices

    virtual int addVertex ( Vertex * vert ) = 0;
    Vertex & vertex ( unsigned int vidx ) const;
    unsigned int nVertices () const;

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
    
    virtual double sourcePhiCoeff ( const double &, const double &,
                                    const double &, const bool,
                                    const std::string & ) const = 0;
    virtual Eigen::Vector3d sourceVCoeff ( const double &, const double &,
                                           const double &, const bool,
                                           const std::string & ) const = 0;
    virtual double doubletPhiCoeff ( const double &, const double &,
                                     const double &, const bool,
                                     const std::string & ) const = 0;
    virtual Eigen::Vector3d doubletVCoeff ( const double &, const double &,
                                            const double &, const bool,
                                            const std::string & ) const = 0;
    virtual double inducedPotential ( const double &, const double &,
                                      const double &, const bool,
                                      const std::string & ) const = 0;
    virtual Eigen::Vector3d inducedVelocity ( const double &, const double &,
                                              const double &, const bool,
                                              const std::string & ) const = 0;

};

#endif
