// Header for node class

#ifndef NODE_H
#define NODE_H

#include <vector>

// Forward declarations

class Edge;
class Face;

/******************************************************************************/
//
// Node class. Defines node location and index, and provides access to connected
// edges and faces.
//
/******************************************************************************/
class Node {

  private:

    unsigned int _lbl;             // Node label
    double _x, _y, _z;             // Vertex coordinates

    unsigned int _curredges, _currfaces;
                                   // # of edges and faces that have been set
                                   
    std::vector<Edge*> _edgeref;   // Vector of edge pointers
    std::vector<Face*> _faceref;   // Vector of face pointers

  public:

    // Constructors

    Node ();
    Node ( unsigned int, const double &, const double &, const double & );
                            // Version with index and coordinates specified

    // Setting and accessing node index

    void setLabel ( unsigned int );
    unsigned int label () const;

    // Setting and accessing vertex coordinates

    void setCoordinates ( const double &, const double &, const double & );
    const double & x () const;
    const double & y () const;
    const double & z () const;

    // Access to edges

    unsigned int numEdges () const;
    void increaseNumEdges ();
    void setNextEdge ( Edge * );
    Edge & edge ( unsigned int ) const;

    // Access to faces

    unsigned int numFaces () const;
    void increaseNumFaces ();
    void setNextFace ( Face * );
    Face & face ( unsigned int ) const;

};

#endif
