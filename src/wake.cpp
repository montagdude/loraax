#include <vector>
#include <Eigen/Core>
#include "util.h"
#include "node.h"
#include "wakeline.h"
#include "wake.h"

/******************************************************************************/
//
// Wake class. Stores wake lines and horseshoe vortices and facilitates their
// setup and connections with each other and the trailing edges of the aircraft
//
/******************************************************************************/

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
Wake::Wake ()
{
  _numlines = 0;
  _numhorseshoes = 0;
}

/******************************************************************************/
//
// Returns number of wake lines
//
/******************************************************************************/
unsigned int Wake::numWakeLines () const { return _numlines; }

/******************************************************************************/
//
// Creates a wake line with given length and number of filaments, aligned with
// the freestream
//
/******************************************************************************/
void Wake::addWakeLine ( const double & length, unsigned int numfilaments,
                         const Eigen::Vector3d & uinf, Node * tenode )
{
  unsigned int i;
  double flen;
  Eigen::Vector3d point;
  WakeLine wakeline;

  // Set number of filaments

  wakeline.setNumFilaments(numfilaments);

  // Set trailing edge node

  wakeline.setTENode(tenode);

  // Compute filament length

  flen = length / double(numfilaments);

  // Set endpoints

  point(0) = tenode->x();
  point(1) = tenode->y();
  point(2) = tenode->z();
  for ( i = 1; i <= numfilaments; i++ )
  {
    point += flen*uinf;
    wakeline.setEndpoint(i, point(0), point(1), point(2));
  }

  // Append to vector of wake lines and update counter

  _wakelines.push_back(wakeline);
  _numlines += 1;
}

/******************************************************************************/
//
// Returns a reference to a wake line
//
/******************************************************************************/
WakeLine & Wake::wakeLine ( unsigned int wlidx )
{
#ifndef NDEBUG

  if (wlidx >= _numlines)
  {
    conditional_stop(1, "Wake::wakeLine", "Invalid wake line specified.");
  }

#endif

  return _wakelines[wlidx];
}
