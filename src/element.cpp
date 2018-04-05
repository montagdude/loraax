#include <vector>
#include "util.h"
#include "vertex.h"
#include "element.h"

/******************************************************************************/
//
// Element class. Can be a panel or vortex element.
// 
/******************************************************************************/

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
Element::Element ()
{
  _idx = -1;
  _currverts = 0;
  _verts.resize(0);
  _type = "";
} 

/******************************************************************************/
//
// Set/access index
//
/******************************************************************************/
void Element::setIdx ( int idx ) { _idx = idx; }
int Element::idx () const { return _idx; }

/******************************************************************************/
//
// Set/access type
//
/******************************************************************************/
void Element::setType ( const std::string & type ) { _type = type; }
const std::string & Element::type () const { return _type; }

/******************************************************************************/
//
// Access vertices
//
/******************************************************************************/
Vertex & Element::vertex ( unsigned int vidx ) const
{
  if (vidx >= _currverts)
  {
    conditional_stop(1, "Element::vertex", "Index out of range.");
  }

  return *_verts[vidx];
}
unsigned int Element::nVertices () const { return _currverts; }
