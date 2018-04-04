// Header for wing class

#ifndef WING_H
#define WING_H

#include <vector>
#include <string>
#include "section.h"
#include "airfoil.h"
#include "quadface.h"
#include "triface.h"
#include "wake.h"

class Vertex;

/******************************************************************************/
//
// Wing class. Contains sections, airfoils, paneling info, faces, a wake, etc.
//
/******************************************************************************/
class Wing {

  private:

    std::string _name;
    unsigned int _nchord, _nspan;	// Number of points along chord and
					//   span. nchord is on each side (top
					//   and bottom); total is 2*nchord-1.
    double _lesprat, _tesprat;		// LE and TE spacing ratios relative to
					//   uniform
    double _rootsprat, _tipsprat;	// Root and tip spacing ratios relative
					//   to uniform

    std::vector<Section> _sections;	// Spanwise sections used for
					//   calculations, set by _nspan and
					//   root & tip spacing rations
    std::vector<Airfoil> _foils;  	// User-specified airfoils
    std::vector<Vertex *> _verts;	// Pointers to vertices on wing and wake
    std::vector<QuadFace> _quads;	// Quad faces (panels)
    std::vector<TriFace> _tris;		// Tri faces (panels) at tip 
    Wake _wake;				// Wake

    std::vector<double> adjustSpacing ( 
                               const std::vector<double> & nom_stations ) const;

  public:

    // Constructor

    Wing ();

    // Set/get name

    void setName ( const std::string & name );
    const std::string & name () const;

    // Set discretization and spacing info

    void setDiscretization ( unsigned int nchord, unsigned int nspan,
                             const double & lesprat, const double & tesprat,
                             const double & rootsprat,
                             const double & tipsprat );

    // Set airfoils. Section airfoils are interpolated from these airfoils.
    // Airfoils do not need to be given in sorted order.

    void setAirfoils ( std::vector<Airfoil> & foils );

    // Set sections based on user section inputs and spacing. User sections
    // define the planform shape, and the _sections variables defines the final
    // discretization. Sections do not need to be given in sorted order.

    int setupSections ( std::vector<Section> & user_sections );

    // Creates panels and surface vertex pointers

    void createPanels ( int & next_global_vertidx, int & next_global_faceid );

    // Set up wake

    void setupWake ( int & next_wake_vertidx, int & next_wake_ringidx );

    // Access to verts and panels

    unsigned int nVerts () const;
    unsigned int nQuads () const;
    unsigned int nTris () const;
    Vertex * vert ( unsigned int vidx );
    QuadFace * quadFace ( unsigned int qidx );
    TriFace * triFace ( unsigned int tidx );

    // Access to wake

    Wake & wake ();
};

#endif
