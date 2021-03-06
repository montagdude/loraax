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

    // Note: all geometric quantities use incompressible (Prandtl-Glauert
    // transformed) coordinates unless suffixed with 'comp'

    double _sigma, _mu;             // Source and doublet strength
    double _length;                 // A characteristic length
    double _area;                   // Face area
    Eigen::Vector3d _norm;          // Outward-facing normal
    Eigen::Vector3d _cen;           // Face centroid
    Eigen::Vector3d _colloc;        // Collocation point (defaults to centroid)
    bool _colloc_is_centroid;       // Whether collocation point is centroid
    Eigen::Matrix3d _trans, _invtrans;
                                    // Transform from inertial frame to panel
                                    // frame and vice versa
    std::vector<double> _xtrans, _ytrans;
                                    // Panel endpoint coordinates in panel frame

    // Geometric quantities in actual (non P-G transformed) coordinates

    double _areacomp;
    Eigen::Vector3d _normcomp, _cencomp;
    Eigen::Vector3d _tancomp;       // Surface tangent vector in streamwise dir.

    Eigen::Vector3d _vel, _velcomp; // Incompressible and compressible velocity
                                    //   at centroid
    double _mach;                   // Local mach number
    double _cp, _p;                 // Pressure coefficient and pressure
    double _rho;                    // Density
    double _cf, _mdefect, _dmdefect;// Skin friction coefficient, mass defect,
                                    //   and d/ds(mass defect)

    Panel * _right, * _left, * _front, * _back;
                                    // Neighbor panels
    Eigen::Matrix3d _jac;           // Grid metrics jacobian matrix
    Eigen::PartialPivLU<Eigen::Matrix3d> _lu;
                                    // Factorization of jacobian matrix 
    
    const static double _farfield_distance_factor;
    
    // Computing geometric quantities. Note: computes both actual and
    // incompressible (Prandtl-Glauert) quantities.
    
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
    
    // Computing and accessing source strength

    void setSourceStrength ( const double & sigma );
    void computeSourceStrength ( const Eigen::Vector3d & uinfvec,
                                 bool viscous );
    const double & sourceStrength () const;

    // Setting and accessing doublet strength

    void setDoubletStrength ( const double & mu );
    const double & doubletStrength () const;

    // Returns geometric quantities. Compressible and P-G coordinates as
    // inidicated

    const double & area () const;
    const double & areaComp () const;
    const Eigen::Vector3d & centroid () const;
    const Eigen::Vector3d & centroidComp () const;
    const Eigen::Vector3d & normal () const;
    const Eigen::Vector3d & normalComp () const;
    const Eigen::Vector3d & tanComp () const;

    // Set surface tangent direction

    void setTangentComp ( const Eigen::Vector3d & tancomp );

    // Set/access collocation point. If not set, defaults to centroid.

    void setCollocationPoint ( const Eigen::Vector3d & colloc );
    const Eigen::Vector3d & collocationPoint () const;
    bool collocationPointIsCentroid () const;

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
    
    // Returns induced velocity due to modeling the doublet panel as a vortex
    // ring, which should give the same result as the doublet panel, but is
    // better behaved and cheaper to compute. Does not include the source
    // contribution.
    
    virtual Eigen::Vector3d vortexVelocity ( const double & x,
                                             const double & y, const double & z,
                                             const double & rcore,
                                             bool mirror_y=false ) const = 0;
    
    // Compute or access surface velocity. Must set neighbors and compute grid
    // transformation before computing velocity.
    
    void computeVelocity ( const Eigen::Vector3d & uinfvec );
    const Eigen::Vector3d & velocity () const;
    const Eigen::Vector3d & velocityComp () const;

    // Mach number

    const double & mach () const;
    
    // Compute or access pressure, density, and pressure coefficient
    
    int computePressure ( const double & uinf, const double & rhoinf,
                          const double & pinf );
    const double & pressure () const;
    const double & density () const;
    const double & pressureCoefficient () const;
    
    // Mass defect: uedge*deltastar and derivative d/ds

    const double & massDefect() const;
    const double & massDefectDerivative () const;

    // Average viscous quantities from vertices

    void averageFromVertices ();

    // Compute force and moment contributions
    
    void computeForceMoment ( const double & uinf, const double & rhoinf,
                              const double & pinf,
                              const Eigen::Vector3d & moment_center,
                              bool viscous, Eigen::Vector3d & fp,
                              Eigen::Vector3d & fv, Eigen::Vector3d & mp,
                              Eigen::Vector3d & mv );
};

#endif
