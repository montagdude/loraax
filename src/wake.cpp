#define _USE_MATH_DEFINES

#include <vector>
#include <cmath>
#include "util.h"
#include "settings.h"
#include "vertex.h"
#include "vortex_ring.h"
#include "wake.h"

/******************************************************************************/
//
// Wake class.
//
/******************************************************************************/

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
Wake::Wake ()
{
  _verts.resize(0);
  _vrings.resize(0); 
}

/******************************************************************************/
//
// initialize with a vector of vertices along the wing trailing edge
//
/******************************************************************************/
void Wake::initialize ( const std::vector<Vertex> & teverts,
                        int & next_wake_vertidx, int & next_wake_ringidx )
{
  int nstream, nspan, nvring;
  unsigned int i, j, blvert, brvert, trvert, tlvert;
  std::vector<double> uinfdir;
  double x, y, z;

  // Determine number of rows to store based on inputs

  nspan = teverts.size();
  nstream = int(ceil(wakelen/(dt*uinf)));
  _verts.resize(nspan*nstream);

  // Add vertices along freestream direction

  uinfdir.resize(3);
  uinfdir[0] = cos(alpha*M_PI/180.);
  uinfdir[1] = 0.;
  uinfdir[2] = sin(alpha*M_PI/180.);

  for ( i = 0; int(i) < nspan; i++ )
  {
    for ( j = 0; int(j) < nstream; j++ )
    {
      if (j == 0)
        _verts[i*nstream] = teverts[i];
      else   
      {
        x = _verts[i*nstream].x() + uinfdir[0]*double(j)*dt*uinf;
        y = _verts[i*nstream].y();
        z = _verts[i*nstream].z() + uinfdir[2]*double(j)*dt*uinf;
        _verts[i*nstream+j].setCoordinates(x, y, z);
      }
      _verts[i*nstream+j].setIdx(next_wake_vertidx);
      next_wake_vertidx += 1;
    }
  }

  // Create vortex rings

  nvring = (nspan-1)*(nstream-1);
  _vrings.resize(nvring);
  for ( i = 0; int(i) < nspan-1; i++ )
  {
    for ( j = 0; int(j) < nstream-1; j++ )
    {
      blvert = i*nstream+j+1;
      brvert = (i+1)*nstream+j+1;
      trvert = (i+1)*nstream+j;
      tlvert = i*nstream+j;
      _vrings[i*(nstream-1)+j].setIdx(next_wake_ringidx);
      _vrings[i*(nstream-1)+j].addVertex(&_verts[blvert]);
      _vrings[i*(nstream-1)+j].addVertex(&_verts[brvert]);
      _vrings[i*(nstream-1)+j].addVertex(&_verts[trvert]);
      _vrings[i*(nstream-1)+j].addVertex(&_verts[tlvert]);
      next_wake_ringidx += 1;
    }
  }
}

/******************************************************************************/
//
// Access to verts, vortex rings
//
/******************************************************************************/
unsigned int Wake::nVerts () const { return _verts.size(); }
unsigned int Wake::nVRings () const { return _vrings.size(); }
Vertex * Wake::vert ( unsigned int vidx )
{
#ifdef DEBUG
  if (vidx >= _verts.size())
    conditional_stop(1, "Wake::vert", "Index out of range.");
#endif

  return &_verts[vidx];
}

VortexRing * Wake::vRing ( unsigned int vridx )
{
#ifdef DEBUG
  if (vridx >= _vrings.size())
    conditional_stop(1, "Wake::vRing", "Index out of range.");
#endif

  return &_vrings[vridx];
}
