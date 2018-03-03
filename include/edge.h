// Header for Face class

#ifndef EDGE_H
#define EDGE_H

#include <vector>

// Forward declarations

class Node;
class Face;

/******************************************************************************/
//
// Edge class.  Stores references to connected faces and nodes, etc.
//
/******************************************************************************/
class Edge {

  private:

    unsigned int _iedge;         // Edge index
    bool _teflag;                // Flag for lifting surface trailing edge

    unsigned int _currfaces;     // # of faces that have been set

    std::vector<Node*> _noderef; // Vector of node pointers (2 per edge) 
    std::vector<Face*> _faceref; // Vector of face pointers (2 per edge)

  public:

    // Constructor

    Edge ();

    // Setting and accessing edge index

    void setIndex ( unsigned int );
    unsigned int index () const;

    // Access to nodes

    void setNode ( unsigned int, Node * );
    Node & node ( unsigned int ) const;

    // Access to faces

    unsigned int numFaces () const;
    void increaseNumFaces ();
    void setNextFace ( Face * );
    Face & face ( unsigned int ) const;

    // Setting and querying whether edge belongs to lifting surface TE

    void setTrailingEdge ();
    bool isTrailingEdge () const;

};

#endif
