#include <vector>
#include "util.h"
#include "node.h"
#include "face.h"
#include "edge.h"

/******************************************************************************/
//
// Edge class.  Stores references to connected faces and nodes, etc.
//
/******************************************************************************/

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
Edge::Edge ()
{
  _teflag = false;
  _noderef.resize(2);
  _currfaces = 0;
}

/******************************************************************************/
//
// Sets edge index
//
/******************************************************************************/
void Edge::setIndex ( unsigned int idx ) { _iedge = idx; }

/******************************************************************************/
//
// Returns edge index
//
/******************************************************************************/
unsigned int Edge::index () const { return _iedge; }

/******************************************************************************/
//
// Sets a pointer to a connected node
//
/******************************************************************************/
void Edge::setNode ( unsigned int nidx, Node * nodein )
{
#ifndef NDEBUG

  if (nidx >= 2)
  {
    conditional_stop(1, "Edge::setNode", "Index must be 0 or 1.");
  }

#endif

  _noderef[nidx] = nodein;
}

/******************************************************************************/
//
// Returns node reference
//
/******************************************************************************/
Node & Edge::node ( unsigned int nidx ) const
{
#ifndef NDEBUG

  if (nidx >= 2)
  {
    conditional_stop(1, "Edge::node", "Index must be 0 or 1.");
  }

#endif

  return *_noderef[nidx];
}

/******************************************************************************/
//
// Returns number of faces currently set
//
/******************************************************************************/
unsigned int Edge::numFaces () const { return _currfaces; }

/******************************************************************************/
//
// Increases size of face vector
//
/******************************************************************************/
void Edge::increaseNumFaces () 
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
void Edge::setNextFace ( Face * facein )
{
#ifndef NDEBUG

  if (_currfaces >= _faceref.size())
  {
    conditional_stop(1, "Edge::setNextFace", 
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
Face & Edge::face ( unsigned int fidx ) const
{
#ifndef NDEBUG

  if (fidx >= _currfaces)
  {
    conditional_stop(1, "Edge::face", "Face index out of bounds.");
  }

#endif

  return *_faceref[fidx];
}

/******************************************************************************/
//
// Sets edge as a lifting surface trailing edge (from which vortex wake will be
// shed)
//
/******************************************************************************/
void Edge::setTrailingEdge () { _teflag = true; }

/******************************************************************************/
//
// Notifies whether edge is a lifting surface trailing edge
//
/******************************************************************************/
bool Edge::isTrailingEdge () const { return _teflag; }
