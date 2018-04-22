#define _USE_MATH_DEFINES

#include <vector>
#include <cmath>
#include "util.h"
#include "settings.h"
#include "vertex.h"
#include "tripanel.h"
#include "quadpanel.h"
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
  _tris.resize(0); 
  _quads.resize(0); 
}

/******************************************************************************/
//
// initialize with a vector of vertices along the wing trailing edge
//
/******************************************************************************/
void Wake::initialize ( const std::vector<Vertex> & teverts,
                        int & next_global_vertidx, int & next_global_elemidx )
{
  int nstream, nspan, ntris, nquads;
  unsigned int i, j, blvert, brvert, trvert, tlvert;
  double x, y, z, xviz, yviz, zviz;

  // Determine number of rows to store based on inputs

  nspan = teverts.size();
  nstream = int(ceil(rollupdist/(dt*uinf))); // Number of relaxable streamwise
                                             // points. One more will be added
                                             // for the trailing horseshoe.
  _verts.resize(nspan*(nstream+1));

  // Add vertices along freestream direction
//FIXME: allow user-defined initial wake angle, because if it is too high,
// solution can go unstable. Checked one case with full pivoting LU instead of
// partial pivoting, but it was still unstable.

  for ( i = 0; int(i) < nspan; i++ )
  {
    for ( j = 0; int(j) < nstream+1; j++ )
    {
      if (j == 0)
        _verts[i*(nstream+1)] = teverts[i];
      else if (int(j) < nstream)
      {
        // Vertices that will roll up
        x = _verts[i*(nstream+1)].x() + uinfvec(0)*double(j)*dt
        y = _verts[i*(nstream+1)].y();
        z = _verts[i*(nstream+1)].z() + uinfvec(2)*double(j)*dt
        _verts[i*(nstream+1)+j].setCoordinates(x, y, z);
      }
      else
      {
        // Trailing vertex at "infinity"
        x = _verts[i*(nstream+1)].x() + uinfvec(0)*1001.*rollupdist/uinf;
        y = _verts[i*(nstream+1)].y();
        z = _verts[i*(nstream+1)].z() + uinfvec(2)*1001.*rollupdist/uinf;
        _verts[i*(nstream+1)+j].setCoordinates(x, y, z);
        xviz = _verts[i*(nstream+1)].x() + uinfvec(0)*nstream*dt;
        yviz = _verts[i*(nstream+1)].y();
        zviz = _verts[i*(nstream+1)].z() + uinfvec(2)*nstream*dt;
        _verts[i*(nstream+1)+j].setVizCoordinates(xviz, yviz, zviz);
      }
      _verts[i*(nstream+1)+j].setIdx(next_global_vertidx);
      next_global_vertidx += 1;
    }
  }

  // Create tri doublets (which can roll up). Tri doublets are used for this
  // portion of the wake because doublet panels are required to be planar, and
  // with quads there would be significant distortion.

  ntris = (nspan-1)*(nstream-1)*2;
  _tris.resize(ntris);
  for ( i = 0; int(i) < nspan-1; i++ )
  {
    for ( j = 0; int(j) < nstream-1; j++ )
    {
      blvert = i*(nstream+1)+j+1;
      brvert = (i+1)*(nstream+1)+j+1;
      trvert = (i+1)*(nstream+1)+j;
      tlvert = i*(nstream+1)+j;

      _tris[i*(nstream-1)*2+j*2].setIdx(next_global_elemidx);
      _tris[i*(nstream-1)*2+j*2].addVertex(&_verts[tlvert]);
      _tris[i*(nstream-1)*2+j*2].addVertex(&_verts[blvert]);
      _tris[i*(nstream-1)*2+j*2].addVertex(&_verts[brvert]);
      next_global_elemidx += 1;

      _tris[i*(nstream-1)*2+j*2+1].setIdx(next_global_elemidx);
      _tris[i*(nstream-1)*2+j*2+1].addVertex(&_verts[brvert]);
      _tris[i*(nstream-1)*2+j*2+1].addVertex(&_verts[trvert]);
      _tris[i*(nstream-1)*2+j*2+1].addVertex(&_verts[tlvert]);
      next_global_elemidx += 1;
    }
  }

  // Create quad doublets which extend to 1000*rollupdist

  nquads = nspan-1;
  _quads.resize(nquads);
  for ( i = 0; int(i) < nspan-1; i++ )
  {
    blvert = i*(nstream+1)+nstream;
    tlvert = i*(nstream+1)+nstream-1;
    trvert = (i+1)*(nstream+1)+nstream-1;
    brvert = (i+1)*(nstream+1)+nstream;
    _quads[i].setIdx(next_global_elemidx);
    _quads[i].addVertex(&_verts[tlvert]);
    _quads[i].addVertex(&_verts[blvert]);
    _quads[i].addVertex(&_verts[brvert]);
    _quads[i].addVertex(&_verts[trvert]);
    next_global_elemidx += 1;
  }
}

/******************************************************************************/
//
// Access to verts, vortex rings
//
/******************************************************************************/
unsigned int Wake::nVerts () const { return _verts.size(); }
unsigned int Wake::nTris () const { return _tris.size(); }
unsigned int Wake::nQuads () const { return _quads.size(); }
Vertex * Wake::vert ( unsigned int vidx )
{
#ifdef DEBUG
  if (vidx >= _verts.size())
    conditional_stop(1, "Wake::vert", "Index out of range.");
#endif

  return &_verts[vidx];
}

TriPanel * Wake::triPanel ( unsigned int tidx )
{
#ifdef DEBUG
  if (tidx >= _tris.size())
    conditional_stop(1, "Wake::triPanel", "Index out of range.");
#endif

  return &_tris[tidx];
}

QuadPanel * Wake::quadPanel ( unsigned int qidx )
{
#ifdef DEBUG
  if (qidx >= _quads.size())
    conditional_stop(1, "Wake::quadPanel", "Index out of range.");
#endif

  return &_quads[qidx];
}
