#include <vector>
#include <Eigen/Core>
#include "util.h"
#include "vertex.h"
#include "singularities.h"
#include "vortex_ring.h"

/******************************************************************************/
//
// Vortex ring class
//
/******************************************************************************/

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
VortexRing::VortexRing ()
{
  _idx = -1;
  _currverts = 0;
  _gamma = 0.;
  _verts.resize(4);
}

/******************************************************************************/
//
// Set/access index
//
/******************************************************************************/
void VortexRing::setIdx ( int idx ) { _idx = idx; }
int VortexRing::idx () const { return _idx; }

/******************************************************************************/
//
// Set or access vertices
//
/******************************************************************************/
int VortexRing::addVertex ( Vertex * vert )
{
  if (_currverts == 4)
  {
    conditional_stop(1, "VortexRing::addVertex",
                     "4 vertices have already been set.");
    return 1;
  }

  _verts[_currverts] = vert;
  _currverts += 1;

  return 0;
}

Vertex & VortexRing::vertex ( unsigned int vidx ) const
{
  if (vidx >= _currverts)
  {
    conditional_stop(1, "VortexRing::vertex", "Index out of range.");
  }

  return *_verts[vidx];
}

/******************************************************************************/
//
// Set or access circulation
//
/******************************************************************************/
void VortexRing::setCirculation ( const double & gamma ) { _gamma = gamma; }
const double & VortexRing::circulation () const { return _gamma; }

/******************************************************************************/
//
// Computes induced velocity at a point due to wake line with unit circulation
//
/******************************************************************************/
Eigen::Vector3d VortexRing::VCoeff ( const double & x, const double & y, 
                                     const double & z,
                                     const double & rcore ) const
{
  unsigned int i;
  Eigen::Vector3d vel;

#ifdef DEBUG
  if (_currverts < 4)
    conditional_stop(1, "VortexRing::VCoeff",
                     "4 vertices have not yet been set.");
#endif

  // Initialize induced velocity

  vel(0) = 0.0;
  vel(1) = 0.0;
  vel(2) = 0.0;

  // First _numfilaments-1 treated as finite filaments

  for ( i = 0; i < 3; i++ )
  {
    vel += vortex_velocity(x, y, z, _verts[i]->x(), _verts[i]->y(),
                           _verts[i]->z(), _verts[i+1]->x(), _verts[i+1]->y(),
                           _verts[i+1]->z(), rcore, false);
  }
  vel += vortex_velocity(x, y, z, _verts[3]->x(), _verts[3]->y(),
                         _verts[3]->z(), _verts[0]->x(), _verts[0]->y(),
                         _verts[0]->z(), rcore, false);

  return vel;
}

/******************************************************************************/
//
// Computes induced velocity at a point
//
/******************************************************************************/
Eigen::Vector3d VortexRing::inducedVelocity ( const double & x,
                                             const double & y, const double & z,
                                             const double & rcore ) const
{
  Eigen::Vector3d vel;

  vel = VCoeff(x, y, z, rcore)*_gamma; 
  return vel;
}
