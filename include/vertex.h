// Header for Vertex class

#ifndef VERTEX_H
#define VERTEX_H

class Face;

/******************************************************************************/
//
// Vertex class. Defines x, y, z coordinates in space and stores references to
// neighboring faces.
//
/******************************************************************************/
class Vertex {

  private:

    double _x, _y, _z;
    std::vector<Face *> _faces;		// Neighboring faces
    unsigned int _nfaces;

  public:

    // Constructor

    Vertex ();

    // Setting or accessing coordinates

    void setCoordinates ( const double & x, const double & y,
                          const double & z );
    const double & x () const;
    const double & y () const;
    const double & z () const;

    // Adding or accessing face references

    int addFace ( Face * face ); 
    Face & face ( unsigned int fidx ) const;
    bool isNeighbor ( const Face * face ) const;
    unsigned int nFaces () const;
};

#endif
