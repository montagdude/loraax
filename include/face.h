// Header for face class

#ifndef FACE_H
#define FACE_H

#include <vector>
#include <string>
#include <Eigen/Core>

/******************************************************************************/
//
// Face class. Computes geometric quantities, source/doublet influence
// coefficients at a point, etc.
//
/******************************************************************************/
class Face {

  protected:

    double _sigma, _mu;            // Source and doublet strength
    double _length;                // A characteristic length
    double _area;                  // Face area
    Eigen::Vector3d _norm;         // Outward-facing normal
    Eigen::Vector3d _cen;          // Face centroid
    Eigen::Matrix3d _trans, _invtrans;
                                   // Transform from inertial frame to panel
                                   // frame and vice versa
    std::vector<double> _x, _y, _z;
                                   // Panel endpoint coordinates
    std::vector<double> _xtrans, _ytrans;
                                   // Panel endpoint coordinates in panel frame

  public:

    // Constructor

    Face ();

    // Initialize geometry with set of endpoints

    virtual void setEndpoints ( const std::vector<double> & x,
                                const std::vector<double> & y,
                                const std::vector<double> & z ) = 0;

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
