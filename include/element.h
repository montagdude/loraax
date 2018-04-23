// Header for element class

#ifndef ELEMENT_H
#define ELEMENT_H

#include <vector>
#include <string>

class Vertex;

/******************************************************************************/
//
// Element class. Can be a panel or vortex element.
//
/******************************************************************************/
class Element {

  protected:

    int _idx;		   	   // Element identifier index
    unsigned int _currverts;
    std::vector<Vertex *> _verts;  // Element endpoint vertices
    std::string _type;		   // Element type

  public:

    // Constructor

    Element ();

    // Set / access index

    void setIdx ( int idx );
    int idx () const;

    // Set / access type

    void setType ( const std::string & type );
    const std::string & type () const;

    // Set / access vertices

    virtual int addVertex ( Vertex * vert, bool ref_element_to_vert=true ) = 0;
    Vertex & vertex ( unsigned int vidx ) const;
    unsigned int nVertices () const;
};

#endif
