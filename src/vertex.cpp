#include <vector>
#include <Eigen/Core>
#include "util.h"
#include "element.h"
#include "vertex.h"

/******************************************************************************/
//
// Vertex class. Defines x, y, z coordinates in space and stores references to
// neighboring elements.
// 
/******************************************************************************/

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
Vertex::Vertex ()
{
  _idx = -1;
  _x = 0.0;
  _y = 0.0;
  _z = 0.0;
  _nelems = 0;
  _elements.resize(0);
  _data.resize(8);
}

/******************************************************************************/
//
// Set/access index
//
/******************************************************************************/
void Vertex::setIdx ( int idx ) { _idx = idx; }
int Vertex::idx () const { return _idx; }

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
// Transformations
//
/******************************************************************************/
void Vertex::scale ( const double & factor )
{
  _x *= factor;
  _y *= factor;
  _z *= factor;
}

void Vertex::translate ( const double & dx, const double & dy,
                         const double & dz )
{
  _x += dx;
  _y += dy;
  _z += dz;
}

void Vertex::rotate ( const Eigen::Matrix3d & transform )
{
  Eigen::Vector3d point, transpoint;

  point(0) = _x;
  point(1) = _y;
  point(2) = _z;
  transpoint = transform*point;
  _x = transpoint(0);
  _y = transpoint(1);
  _z = transpoint(2);
}

/******************************************************************************/
//
// Adding or accessing element references
//
/******************************************************************************/
int Vertex::addElement ( Element * element )
{
  unsigned int i;

  // Only add if the element is not already in the list

  for ( i = 0; i < _nelems; i++ )
  {
    if (_elements[i]->idx() == element->idx())
    {
#ifdef DEBUG
      print_warning("Vertex::addElement", "Element " + 
                    int2string(element->idx()) +
                    std::string(" already in list."));
#endif
      return 1;
    }
  }

  _nelems += 1;
  _elements.push_back(element);

  return 0;
}

Element * Vertex::element ( unsigned int eidx )
{
  if (eidx >= _nelems)
  {
    conditional_stop(1, "Vertex::element", "Index out of range.");
  }

  return _elements[eidx];
}

bool Vertex::isNeighbor ( const Element * element ) const
{
  unsigned int i;

  for ( i = 0; i < _nelems; i++ )
  {
    if (_elements[i]->idx() == element->idx())
      return true;
  }

  return false;
}

unsigned int Vertex::nElems () const { return _nelems; }

/******************************************************************************/
//
// Setting or accessing data. See legend in comments in vertex.h.
//
/******************************************************************************/
int Vertex::setData ( unsigned int idx, const double & var )
{
#ifdef DEBUG
  if (idx >= _data.size())
  {
    conditional_stop(1, "Vertex::setData", "Index out of range.");
    return 1;
  }
#endif

  _data[idx] = var;

  return 0;
}

const double & Vertex::data ( unsigned int idx ) const
{
#ifdef DEBUG
  if (idx >= _data.size())
    conditional_stop(1, "Vertex::data", "Index out of range.");
#endif

  return _data[idx];
}
