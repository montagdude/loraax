#include <vector>
#include "util.h"
#include "face.h"
#include "vertex.h"

/******************************************************************************/
//
// Vertex class. Defines x, y, z coordinates in space and stores references to
// neighboring faces.
// 
/******************************************************************************/

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
Vertex::Vertex ()
{
  _x = 0.0;
  _y = 0.0;
  _z = 0.0;
  _nfaces = 0;
  _faces.resize(0);
}

/******************************************************************************/
//
// Set or access coordinates
//
/******************************************************************************/
void Vertex::setCoordinates ( const double & x, const double & y,
                              const double & z )
{
  _x = x;
  _y = y;
  _z = z;
}

const double & Vertex::x () const { return _x; }
const double & Vertex::y () const { return _y; }
const double & Vertex::z () const { return _z; }

/******************************************************************************/
//
// Adding or accessing face references
//
/******************************************************************************/
int Vertex::addFace ( Face * face )
{
  unsigned int i;

  // Only add if the face is not already in the list

  for ( i = 0; i < _nfaces; i++ )
  {
    if (_faces[i]->idx() == face->idx())
    {
#ifndef NDEBUG
      print_warning("Vertex::addFace", "Face " + int2string(face->idx()) +
                    std::string(" already in list."));
#endif
      return 1;
    }
  }

  _nfaces += 1;
  _faces.push_back(face);

  return 0;
}

Face & Vertex::face ( unsigned int fidx ) const
{
  if (fidx >= _nfaces)
  {
    conditional_stop(1, "Vertex::face", "Index out of range.");
  }

  return *_faces[fidx];
}

bool Vertex::isNeighbor ( const Face * face ) const
{
  unsigned int i;

  for ( i = 0; i < _nfaces; i++ )
  {
    if (_faces[i]->idx() == face->idx())
      return true;
  }

  return false;
}

unsigned int Vertex::nFaces () const { return _nfaces; }
