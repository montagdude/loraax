#include <vector>
#include <Eigen/Core>
#include "util.h"
#include "vertex.h"
#include "face.h"

/******************************************************************************/
//
// Face class. Computes geometric quantities, source/doublet influence
// coefficients at a point, etc.
// 
/******************************************************************************/

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
Face::Face ()
{
  _idx = -1;
  _sigma = 0.0;
  _mu = 0.0;
  _length = 0.0;
  _area = 0.0;
  _currverts = 0;
  _verts.resize(0);
  _xtrans.resize(0);
} 

/******************************************************************************/
//
// Set/access index
//
/******************************************************************************/
void Face::setIdx ( int idx ) { _idx = idx; }
int Face::idx () const { return _idx; }

/******************************************************************************/
//
// Access vertices
//
/******************************************************************************/
Vertex & Face::vertex ( unsigned int vidx ) const
{
  if (vidx >= _currverts)
  {
    conditional_stop(1, "Face::vertex", "Index out of range.");
  }

  return *_verts[vidx];
}
unsigned int Face::nVertices () const { return _currverts; }

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
