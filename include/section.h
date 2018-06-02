// Header for section class

#ifndef SECTION_H
#define SECTION_H

#include <vector>
#include <string>
#include "sectional_object.h"
#include "vertex.h"
#include "airfoil.h"

/******************************************************************************/
//
// Section class. Defines a wing section.
//
/******************************************************************************/
class Section: public SectionalObject {

  private:

    double _xle, _zle;	// Leading edge location
    double _chord;
    double _twist;	// Twist angle, positive leading edge up
    double _roll;	// Roll angle, positive CCW viewed from back

    unsigned int _nverts;
    std::vector<Vertex>	_verts, _uverts;
			// Vertices defining panel endpoints, and non-rotated,
                        // non-translated version of same
    Airfoil _foil;	// Airfoil at this section
    double _re;         // Reynolds number
    double _fa, _fn;    // Axial and normal force components in section frame
    double _cl, _cd;    // Sectional lift and drag coefficients
    bool _converged;    // Whether Xfoil BL calculations converged

  public:

    // Constructor

    Section (); 

    // Set or access position, orientation, and scale

    void setGeometry ( const double & xle, const double & y, const double & zle,
                       const double & chord, const double & twist,
                       const double & roll );
    void setRoll ( const double & roll );
    const double & xle () const;
    const double & zle () const;
    const double & chord () const;
    const double & twist () const;
    const double & roll () const;

    // Access vertices

    Vertex & vert ( unsigned int idx );

    // Set vertices from spacing distribution

    void setVertices ( unsigned int nchord, const double & lesprat,
                       const double & tesprat );

    // Access airfoil

    Airfoil & airfoil ();

    // Computes section pressure force

    void computePressureForce ( const double & alpha, const double & uinf,
                                const double & rhoinf );

    // BL calculations with Xfoil

    void computeBL ( const Eigen::Vector3d & uinfvec, const double & rhoinf ); 
    bool blConverged () const;

    // Sectional lift and drag coefficients

    const double & liftCoefficient () const;
    const double & dragCoefficient () const;
    
    // Computes Reynolds number and sets it for airfoil
    
    void computeReynoldsNumber ( const double & rhoinf, const double & uinf,
                                 const double & muinf );
    const double & reynoldsNumber () const; 
    
    // Sets Mach number for airfoil
    
    void setMachNumber ( const double & mach );
};

#endif
