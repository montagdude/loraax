#include <vector>
#include <string>
#include <Eigen/Core>
#include "util.h"
#include "node.h"
#include "edge.h"
#include "face.h"

/******************************************************************************/
//
// Face class.  Stores references to connected edges and nodes, computes
// source/doublet influence coefficients at a point, etc.
// 
/******************************************************************************/

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
Face::Face ()
{
  _lbl = 0;
  _thinflag = false;
  _TEuflag = false;
  _TElflag = false;
  _sigma = 0.0;
  _mu = 0.0;
  _area = 0.0;
  _cen.resize(3);
  _norm.resize(3);
  _trans.resize(3,3);
  _invtrans.resize(3,3);
  _currnodes = 0;
  _curredges = 0;
} 

/******************************************************************************/
//
// Sets face label
//
/******************************************************************************/
void Face::setLabel ( unsigned int lbl ) { _lbl = lbl; }

/******************************************************************************/
//
// Returns face label
//
/******************************************************************************/
unsigned int Face::label () const { return _lbl; }

/******************************************************************************/
//
// Returns number of nodes
//
/******************************************************************************/
unsigned int Face::numNodes () const { return _currnodes; }

/******************************************************************************/
//
// Sets a pointer to a connected node
//
/******************************************************************************/
void Face::setNode ( unsigned int nidx, Node * nodein ) 
{ 
  _noderef[nidx] = nodein;
  _currnodes += 1;
}

/******************************************************************************/
//
// Returns node reference
//
/******************************************************************************/
Node & Face::node ( unsigned int nidx ) const
{
#ifndef NDEBUG

  if (nidx >= _currnodes)
  {
    conditional_stop(1, "Face::node", "Invalid node referenced.");
  }

#endif

  return *_noderef[nidx];
}

/******************************************************************************/
//
// Returns number of edges currently set
//
/******************************************************************************/
unsigned int Face::numEdges () const { return _curredges; }

/******************************************************************************/
//
// Sets a pointer to a connected edge
//
/******************************************************************************/
void Face::setNextEdge ( Edge * edgein ) 
{ 
#ifndef NDEBUG

  if (_curredges >= _edgeref.size())
  {
    conditional_stop(1, "Face::setNextEdge", "Too many edges specified.");
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
Edge & Face::edge ( unsigned int eidx ) const
{
#ifndef NDEBUG

  if (eidx >= _curredges)
  {
    conditional_stop(1, "Face::edge", "Invalid edge referenced.");
  }

#endif

  return *_edgeref[eidx];
}

/******************************************************************************/
//
// Sets source strength
//
/******************************************************************************/
void Face::setSourceStrength ( const double & sigin ) { _sigma = sigin; }

/******************************************************************************/
//
// Returns source strength
//
/******************************************************************************/
const double & Face::sourceStrength () const { return _sigma; }

/******************************************************************************/
//
// Sets doublet strength
//
/******************************************************************************/
void Face::setDoubletStrength ( const double & muin ) { _mu = muin; }

/******************************************************************************/
//
// Returns doublet strength
//
/******************************************************************************/
const double & Face::doubletStrength () const { return _mu; }

/******************************************************************************/
//
// Returns face area
//
/******************************************************************************/
const double & Face::area () const { return _area; }

/******************************************************************************/
//
// Returns centroid
//
/******************************************************************************/
const Eigen::Vector3d & Face::centroid () const { return _cen; }

/******************************************************************************/
//
// Returns normal vector
//
/******************************************************************************/
const Eigen::Vector3d & Face::normal () const { return _norm; }

/******************************************************************************/
//
// Generic querying function for scalars
//
/******************************************************************************/
const double & Face::getScalar ( const std::string & varname ) const
{
  if (varname == "area") { return _area; }
  else 
  { 
    conditional_stop(1, "Face::getScalar", 
                     "Unrecognized output variable " + varname + ".");
  }

  return _area;   // Satisfy compiler warning (never reaches this statement)
}

/******************************************************************************/
//
// Generic querying function for vectors
//
/******************************************************************************/
const Eigen::Vector3d & Face::getVector ( const std::string & varname ) const
{
  if (varname == "normal") { return _norm; }
  else 
  { 
    conditional_stop(1, "Face::getVector", 
                     "Unrecognized output variable " + varname + ".");
  }

  return _norm;   // Satisfy compiler warning (never reaches this statement)
}

/******************************************************************************/
//
// Sets face as part of a thin surface (zero thickness)
//
/******************************************************************************/
void Face::setThin () { _thinflag = true; }

/******************************************************************************/
//
// Notifies whether face is part of a thin surface (zero thickness)
//
/******************************************************************************/
bool Face::isThin () const { return _thinflag; }

/******************************************************************************/
//
// Sets face as being an upper trailing edge element (from which a vortex wake
// will be shed)
//
/******************************************************************************/
void Face::setUpperTrailingEdge ()
{
#ifndef NDEBUG

  if (_TElflag)
  {
    print_warning("Face::setUpperTrailingEdge",
                  "Face is already set as lower trailing edge! Switching.");
    _TElflag = false;
  }

#endif

  _TEuflag = true;
}

/******************************************************************************/
//
// Sets face as being a lower trailing edge element (from which a vortex wake
// will be shed)
//
/******************************************************************************/
void Face::setLowerTrailingEdge ()
{
#ifndef NDEBUG

  if (_TEuflag)
  {
    print_warning("Face::setLowerTrailingEdge",
                  "Face is already set as upper trailing edge! Switching.");
    _TEuflag = false;
  }

#endif

  _TElflag = true;
}

/******************************************************************************/
//
// Notifies whether this is a trailing edge face
//
/******************************************************************************/
bool Face::isTrailingEdge () const
{ 
  if (_TEuflag || _TElflag ) { return true; }
  else { return false; }
}

/******************************************************************************/
//
// Returns type of trailing edge face ("upper", "lower", or "none")
//
/******************************************************************************/
std::string Face::trailingEdgeType () const
{
  std::string tetype;

  if (_TEuflag) { tetype = "upper"; }
  else if (_TElflag) { tetype = "lower"; }
  else { tetype = "none"; }

  return tetype;
}
