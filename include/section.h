// Header for section class

#ifndef SECTION_H
#define SECTION_H

#include <vector>
#include <string>
#include "vertex.h"
#include "airfoil.h"

/******************************************************************************/
//
// Section class. Defines a wing section.
//
/******************************************************************************/
class Section {

  private:

    double _y;
    double _xle, _zle;	// Leading edge location
    double _chord;
    double _twist;	// Twist angle, positive leading edge up
    double _roll;	// Roll angle, positive CCW viewed from back

    unsigned int _nverts;
    std::vector<Vertex>	_verts;
			// Vertices defining panel endpoints
    Airfoil _foil;	// Airfoil at this section

  public:

    // Constructor

    Section (); 

    // Set or access geometry

    void setGeometry ( const double & xle, const double & y, const double & zle,
                       const double & chord, const double & twist,
                       const double & roll );
    const double & y () const;
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
};

#endif
