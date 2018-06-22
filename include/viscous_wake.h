// Header for Viscous Wake class

#ifndef VISCOUSWAKE_H
#define VISCOUSWAKE_H

#include <vector>
#include "tripanel.h"

class Section;
class Vertex;

/******************************************************************************/
//
// Viscous wake class. This is the non-convecting wake which provides the proper
// smooth decrease in displacement thickness behind the trailing edge, which is
// needed for correct pressure drag and smooth pressure distribution near the
// trailing edge in viscous cases. The viscous wake points are set by the
// boundary layer solution.
//
/******************************************************************************/
class ViscousWake {

	private:

	unsigned int _nspan, _nwake;			// # points in span and stream dir.
	std::vector<Vertex *> _verts;			// Pointers to BL wake vertices
	std::vector<TriPanel> _tris;			// Tri source panels

	public:

	// Constructor

	ViscousWake ();

	// Initialize with sections from the wing. Note that the 2D section wake is
	// not initialized until after the BL solution is computed for the first
	// time, so this must be called only after that point.

	void initialize ( std::vector<Section> & sections,
	                  int & next_global_vertidx, int & next_global_elemidx );

	// Updates panels and source strengths

	void update ();

	// Access vertices and panels

	unsigned int nVerts () const;
	unsigned int nTris () const;
	Vertex * vert ( unsigned int vidx );
	TriPanel * triPanel ( unsigned int tidx );
};

#endif
