#include <vector>
#include "util.h"
#include "edge.h"
#include "face.h"
#include "node.h"

/******************************************************************************/
//
// Node class. Defines node location and index, and provides access to connected
// edges and faces.
//
/******************************************************************************/

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
Node::Node ()
{
  _lbl = 0;
  _x = 0.0;
  _y = 0.0;
  _z = 0.0;
  _curredges = 0;
  _currfaces = 0;
}

/******************************************************************************/
//
// Constructor with index and coordinates specified
//
/******************************************************************************/
Node::Node ( unsigned int lbl, const double & x, const double & y, 
             const double & z )
{
  _lbl = lbl;
  _x = x;
  _y = y;
  _z = z;
  _curredges = 0;
  _currfaces = 0;
}

/******************************************************************************/
//
// Sets node label
//
/******************************************************************************/
void Node::setLabel ( unsigned int lbl ) { _lbl = lbl; }

/******************************************************************************/
//
// Returns node label
//
/******************************************************************************/
unsigned int Node::label () const { return _lbl; }

/******************************************************************************/
//
// Sets node coordinates
//
/******************************************************************************/
void Node::setCoordinates ( const double & x, const double & y,
                            const double & z)
{
  _x = x;
  _y = y;
  _z = z;
}

/******************************************************************************/
//
// Returns x coordinate
//
/******************************************************************************/
const double & Node::x () const { return _x; }

/******************************************************************************/
//
// Returns y coordinate
//
/******************************************************************************/
const double & Node::y () const { return _y; }

/******************************************************************************/
//
// Returns z coordinate
//
/******************************************************************************/
const double & Node::z () const { return _z; }

/******************************************************************************/
//
// Returns number of edges
//
/******************************************************************************/
unsigned int Node::numEdges () const { return _edgeref.size(); }

/******************************************************************************/
//
// Increases size of edge vector
//
/******************************************************************************/
void Node::increaseNumEdges () 
{ 
  unsigned int numedges;

  numedges = _edgeref.size();
  _edgeref.resize(numedges + 1);
}

/******************************************************************************/
//
// Sets a pointer to a connected edge
//
/******************************************************************************/
void Node::setNextEdge ( Edge * edgein )
{
#ifndef NDEBUG

  if (_curredges >= _edgeref.size())
  {
    conditional_stop(1, "Node::setNextEdge", 
                     "Too many edges specified. Call increaseNumEdges first.");
  }

#endif

  _edgeref[_curredges] = edgein;
  _curredges += 1;
}

/******************************************************************************/
//
// Returns edge reference
//
/******************************************************************************/
Edge & Node::edge ( unsigned int eidx ) const
{
#ifndef NDEBUG

  if (eidx >= _curredges)
  {
    conditional_stop(1, "Node::edge", "Edge index out of bounds.");
  }

#endif

  return *_edgeref[eidx];
}

/******************************************************************************/
//
// Returns number of faces
//
/******************************************************************************/
unsigned int Node::numFaces () const { return _faceref.size(); }

/******************************************************************************/
//
// Increases size of face vector
//
/******************************************************************************/
void Node::increaseNumFaces () 
{ 
  unsigned int numfaces;

  numfaces = _faceref.size();
  _faceref.resize(numfaces + 1);
}

/******************************************************************************/
//
// Sets a pointer to a connected face
//
/******************************************************************************/
void Node::setNextFace ( Face * facein )
{
#ifndef NDEBUG

  if (_currfaces >= _faceref.size())
  {
    conditional_stop(1, "Node::setNextFace", 
                     "Too many faces specified. Call increaseNumFaces first.");
  }

#endif

  _faceref[_currfaces] = facein;
  _currfaces += 1;
}

/******************************************************************************/
//
// Returns face reference
//
/******************************************************************************/
Face & Node::face ( unsigned int fidx ) const
{
#ifndef NDEBUG

  if (fidx >= _currfaces)
  {
    conditional_stop(1, "Node::face", "Face index out of bounds.");
  }

#endif

  return *_faceref[fidx];
}
