// Header for body class

#ifndef BODY_H
#define BODY_H

#include <vector>
#include "node.h"
#include "edge.h"
#include "face.h"
#include "triface.h"
#include "quadface.h"

/******************************************************************************/
//
// Body class. Stores groups of nodes, edges, and faces making up solid bodies.
//
/******************************************************************************/
class Body {

  private:

    unsigned int _nnodes;         // Number of nodes
    unsigned int _nedges;         // Number of edges
    unsigned int _nfaces;         // Number of faces
    unsigned int _nthin;          // Number of faces belonging to thin surfaces
    unsigned int _nte;            // Number of trailing edges
    unsigned int _ntefaces;       // Number of trailing edge faces
    unsigned int _firsttri, _lasttri, _firstquad, _lastquad;
                                  // _faces array idx of first/last tri/quad
    
    std::vector<Node> _nodes;     // Array of nodes
    std::vector<Edge> _edges;     // Array of nodes
    std::vector<TriFace> _tris;   // Array of tri faces
    std::vector<QuadFace> _quads; // Array of quad faces

    std::vector<Face*> _faces;    // Pointer array of faces (tris and quads)

    // Node index from label

    unsigned int nodeIndexFromLabel ( unsigned int ) const;

    // Edge index from node labels
    
    int edgeFromNodes ( unsigned int, unsigned int ) const;

  public:

    // Constructor

    Body ();

    // Adding and accessing nodes

    void addNode ( unsigned int, const double &, const double &, 
                   const double & );
    Node & node ( unsigned int );

    // Adding and accessing faces (first addFace for tris; second for quads)

    void addFace ( unsigned int, unsigned int, unsigned int, unsigned int,
                   bool );
    void addFace ( unsigned int, unsigned int, unsigned int, unsigned int,
                   unsigned int, bool );
    Face & face ( unsigned int );

    // Sets face pointers after all Tris and Quads are read

    void setFacePointers ();

    // Constructs edges from face and node information

    void constructEdges ();

    // Accessing edges

    Edge & edge ( unsigned int );
      
    // Sets trailing edge from node labels

    void setTrailingEdgeFromNodes ( unsigned int, unsigned int );

    // Number of nodes

    unsigned int numNodes () const;

    // Number of faces

    unsigned int numFaces () const;
    unsigned int numTris () const;
    unsigned int numQuads () const;
    unsigned int numThin () const;
    unsigned int numTEFaces () const;

    // Number of edges

    unsigned int numEdges () const;
    unsigned int numTrailingEdges () const;

    // Loop indices in _faces array

    unsigned int firstTri () const;
    unsigned int lastTri () const;
    unsigned int firstQuad () const;
    unsigned int lastQuad () const;

};

#endif
