// Header for Wake class

#ifndef WAKE_H
#define WAKE_H

#include <vector>
#include "vertex.h"
#include "vortex_ring.h"
#include "horseshoe_vortex.h"

/******************************************************************************/
//
// Wake class. Stores vortex rings and horseshoe vortices.
//
/******************************************************************************/
class Wake {

  private:

    std::vector<Vertex> _verts;			// Vertices
    std::vector<VortexRing> _vrings;		// Vortex rings
    std::vector<HorseshoeVortex> _hshoes;	// Horseshoe vortices

  public:

    // Constructor

    Wake ();

    // Initialize with a vector of vertices along a wing trailing edge

    void initialize ( const std::vector<Vertex> & teverts,
                      int & next_global_vertidx, int & next_global_elemidx );

    // Access vertices, vortex rings, horseshoe vortices

    unsigned int nVerts () const;
    unsigned int nVRings () const;
    unsigned int nHShoes () const;
    Vertex * vert ( unsigned int vidx );
    VortexRing * vRing ( unsigned int vridx );
    HorseshoeVortex * hShoe ( unsigned int hsidx );
};

#endif
