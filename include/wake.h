// Header for Wake class

#ifndef WAKE_H
#define WAKE_H

#include <vector>
#include "vertex.h"
#include "vortex_ring.h"

/******************************************************************************/
//
// Wake class. 
//
/******************************************************************************/
class Wake {

  private:

    std::vector<Vertex> _verts;		// Vertices
    std::vector<VortexRing> _vrings;	// Vortex rings

  public:

    // Constructor

    Wake ();

    // Initialize with a vector of vertices along a wing trailing edge

    void initialize ( const std::vector<Vertex> & teverts,
                      int & next_wake_vertidx, int & next_wake_ringidx );

    // Access vertices, vortex rings

    unsigned int nVerts () const;
    unsigned int nVRings () const;
    Vertex * vert ( unsigned int vidx );
    VortexRing * vRing ( unsigned int vridx );
};

#endif
