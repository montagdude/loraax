// Header for Wake class

#ifndef WAKE_H
#define WAKE_H

#include <vector>
#include "vertex.h"
#include "tripanel.h"
#include "quadpanel.h"

/******************************************************************************/
//
// Wake class. Stores vortex rings and horseshoe vortices.
//
/******************************************************************************/
class Wake {

  private:

    std::vector<Vertex> _verts;			// Vertices
    std::vector<TriPanel> _tris;		// Tri doublet panels
    std::vector<QuadPanel> _quads;		// Quad doublet panels (trailing
                                                //   to near-infinity)

  public:

    // Constructor

    Wake ();

    // Initialize with a vector of vertices along a wing trailing edge

    void initialize ( const std::vector<Vertex> & teverts,
                      int & next_global_vertidx, int & next_global_elemidx );

    // Access vertices and panels

    unsigned int nVerts () const;
    unsigned int nTris () const;
    unsigned int nQuads () const;
    Vertex * vert ( unsigned int vidx );
    TriPanel * triPanel ( unsigned int tidx );
    QuadPanel * quadPanel ( unsigned int qidx );
};

#endif
