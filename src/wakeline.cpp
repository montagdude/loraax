#include <vector>
#include <Eigen/Core>
#include "util.h"
#include "node.h"
#include "singularities.h"
#include "wakeline.h"

/******************************************************************************/
//
// Wake line class.  Represents a collection of wake filament trailers behind
// a trailing edge node.
//
/******************************************************************************/

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
WakeLine::WakeLine ()
{
  _numfilaments = 0;
}

/******************************************************************************/
//
// Sets number of filaments and resizes endpoints vectors
//
/******************************************************************************/
void WakeLine::setNumFilaments ( unsigned int numfilaments )
{
#ifndef NDEBUG

  if (_numfilaments != 0)
  {
    print_warning("WakeLine::setNumFilaments",
                  "Resetting number of filaments.");
  }

#endif

  _numfilaments = numfilaments;
  _x.resize(_numfilaments+1);
  _y.resize(_numfilaments+1);
  _z.resize(_numfilaments+1);
  _xnew.resize(_numfilaments+1);
  _ynew.resize(_numfilaments+1);
  _znew.resize(_numfilaments+1);
}

/******************************************************************************/
//
// Returns number of filaments
//
/******************************************************************************/
unsigned int WakeLine::numFilaments () const { return _numfilaments; }

/******************************************************************************/
//
// Sets endpoint coordinates for a filament
//
/******************************************************************************/
void WakeLine::setEndpoint ( unsigned int idx, const double & x, 
                             const double & y, const double & z )
{
#ifndef NDEBUG

  if (idx == 0)
  {
    conditional_stop(1, "WakeLine::setEndpoint",
                  "First endpoint can only be set via the trailing edge node.");
  }
  else if (idx == _numfilaments+1)
  {
    conditional_stop(1, "WakeLine::setEndpoint", "Invalid endpoint specified.");
  }

#endif

  _x[idx] = x;
  _y[idx] = y;
  _z[idx] = z;
}

/******************************************************************************/
//
// Returns x coordinate of an endpoint
//
/******************************************************************************/
const double & WakeLine::x ( unsigned int idx ) const { return _x[idx]; }

/******************************************************************************/
//
// Returns y coordinate of an endpoint
//
/******************************************************************************/
const double & WakeLine::y ( unsigned int idx ) const { return _y[idx]; }

/******************************************************************************/
//
// Returns z coordinate of an endpoint
//
/******************************************************************************/
const double & WakeLine::z ( unsigned int idx ) const { return _z[idx]; }

/******************************************************************************/
//
// Sets new endpoint coordinates for a filament
//
/******************************************************************************/
void WakeLine::setNewEndpoint ( unsigned int idx, const double & x, 
                                const double & y, const double & z )
{
#ifndef NDEBUG

  if (idx == 0)
  {
    conditional_stop(1, "WakeLine::setNewEndpoint",
                  "First endpoint can only be set via the trailing edge node.");
  }
  else if (idx == _numfilaments+1)
  {
    conditional_stop(1, "WakeLine::setNewEndpoint",
                     "Invalid endpoint specified.");
  }

#endif

  _xnew[idx] = x;
  _ynew[idx] = y;
  _znew[idx] = z;
}

/******************************************************************************/
//
// Updates endpoints by copying new coordinates into old.  Used for relaxing
// the wake.
//
/******************************************************************************/
void WakeLine::updateEndpoints ()
{
  unsigned int i;

  for ( i = 1; i <= _numfilaments; i++ )
  {
    _x[i] = _xnew[i];
    _y[i] = _ynew[i];
    _z[i] = _znew[i];
  }
}

/******************************************************************************/
//
// Sets pointer to trailing edge node. Also sets first endpoint to node 
// coordinates.
//
/******************************************************************************/
void WakeLine::setTENode ( Node * tenode ) 
{ 
  _tenode = tenode; 
  _x[0] = _tenode->x();
  _y[0] = _tenode->y();
  _z[0] = _tenode->z();
  _xnew[0] = _tenode->x();
  _ynew[0] = _tenode->y();
  _znew[0] = _tenode->z();
}

/******************************************************************************/
//
// Returns reference to trailing edge node
//
/******************************************************************************/
Node & WakeLine::TENode () const { return *_tenode; }

/******************************************************************************/
//
// Computes induced velocity at a point due to wake line with unit circulation
//
/******************************************************************************/
Eigen::Vector3d WakeLine::VCoeff ( const double & x, const double & y, 
                                  const double & z, const double & rcore ) const
{
  unsigned int i;
  Eigen::Vector3d vel;

  // Initialize induced velocity

  vel(0) = 0.0;
  vel(1) = 0.0;
  vel(2) = 0.0;

  // First _numfilaments-1 treated as finite filaments

  for ( i = 0; i < _numfilaments-1; i++ )
  {
    vel += vortex_velocity(x, y, z, _x[i], _y[i], _z[i], _x[i+1], _y[i+1],
                           _z[i+1], rcore, false);
  }

  // Last filament treated as semi-infinite

  vel += vortex_velocity(x, y, z, _x[_numfilaments-1], _y[_numfilaments-1],
                         _z[_numfilaments-1], _x[_numfilaments], 
                         _y[_numfilaments], _z[_numfilaments], rcore, true);

  return vel;
}
