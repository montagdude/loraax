#define _USE_MATH_DEFINES

#include <vector>
#include <cmath>
#include "util.h"
#include "settings.h"
#include "vertex.h"
#include "vortex_ring.h"
#include "horseshoe_vortex.h"
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
  _hshoes.resize(0); 
}

/******************************************************************************/
//
// initialize with a vector of vertices along the wing trailing edge
//
/******************************************************************************/
void Wake::initialize ( const std::vector<Vertex> & teverts,
                        int & next_global_vertidx, int & next_global_elemidx )
{
  int nstream, nspan, nvring, nhshoe;
  unsigned int i, j, blvert, brvert, trvert, tlvert;
  std::vector<double> uinfdir;
  double x, y, z;

  // Determine number of rows to store based on inputs

  nspan = teverts.size();
  nstream = int(ceil(wakelen/(dt*uinf))); // Number of relaxable streamwise
                                          // points. One more will be added for
                                          // the trailing horseshoe.
  _verts.resize(nspan*(nstream+1));

  // Add vertices along freestream direction

  uinfdir.resize(3);
  uinfdir[0] = cos(alpha*M_PI/180.);
  uinfdir[1] = 0.;
  uinfdir[2] = sin(alpha*M_PI/180.);

  for ( i = 0; int(i) < nspan; i++ )
  {
    for ( j = 0; int(j) < nstream+1; j++ )
    {
      if (j == 0)
        _verts[i*(nstream+1)] = teverts[i];
      else   
      {
        x = _verts[i*(nstream+1)].x() + uinfdir[0]*double(j)*dt*uinf;
        y = _verts[i*(nstream+1)].y();
        z = _verts[i*(nstream+1)].z() + uinfdir[2]*double(j)*dt*uinf;
        _verts[i*(nstream+1)+j].setCoordinates(x, y, z);
      }
      _verts[i*(nstream+1)+j].setIdx(next_global_vertidx);
      next_global_vertidx += 1;
    }
  }

  // Create vortex rings

  nvring = (nspan-1)*(nstream-1);
  _vrings.resize(nvring);
  for ( i = 0; int(i) < nspan-1; i++ )
  {
    for ( j = 0; int(j) < nstream-1; j++ )
    {
      blvert = i*(nstream+1)+j+1;
      brvert = (i+1)*(nstream+1)+j+1;
      trvert = (i+1)*(nstream+1)+j;
      tlvert = i*(nstream+1)+j;
      _vrings[i*(nstream-1)+j].setIdx(next_global_elemidx);
      _vrings[i*(nstream-1)+j].addVertex(&_verts[blvert]);
      _vrings[i*(nstream-1)+j].addVertex(&_verts[tlvert]);
      _vrings[i*(nstream-1)+j].addVertex(&_verts[trvert]);
      _vrings[i*(nstream-1)+j].addVertex(&_verts[brvert]);
      next_global_elemidx += 1;
    }
  }

  // Create horseshoe vortices

  nhshoe = nspan-1;
  _hshoes.resize(nhshoe);
  for ( i = 0; int(i) < nspan-1; i++ )
  {
    blvert = i*(nstream+1)+nstream;
    tlvert = i*(nstream+1)+nstream-1;
    trvert = (i+1)*(nstream+1)+nstream-1;
    brvert = (i+1)*(nstream+1)+nstream;
    _hshoes[i].setIdx(next_global_elemidx);
    _hshoes[i].addVertex(&_verts[blvert]);
    _hshoes[i].addVertex(&_verts[tlvert]);
    _hshoes[i].addVertex(&_verts[trvert]);
    _hshoes[i].addVertex(&_verts[brvert]);
    next_global_elemidx += 1;
  }
}

/******************************************************************************/
//
// Access to verts, vortex rings
//
/******************************************************************************/
unsigned int Wake::nVerts () const { return _verts.size(); }
unsigned int Wake::nVRings () const { return _vrings.size(); }
unsigned int Wake::nHShoes () const { return _hshoes.size(); }
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

HorseshoeVortex * Wake::hShoe ( unsigned int hsidx )
{
#ifdef DEBUG
  if (hsidx >= _hshoes.size())
    conditional_stop(1, "Wake::hShoe", "Index out of range.");
#endif

  return &_hshoes[hsidx];
}
