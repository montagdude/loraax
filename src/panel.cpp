#include <vector>
#include <Eigen/Core>
#include "util.h"
#include "vertex.h"
#include "element.h"
#include "panel.h"

/******************************************************************************/
//
// Panel class. Computes geometric quantities, source/doublet influence
// coefficients at a point, etc.
// 
/******************************************************************************/

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
Panel::Panel ()
{
  _sigma = 0.0;
  _mu = 0.0;
  _length = 0.0;
  _area = 0.0;
  _xtrans.resize(0);
  _vel << 0., 0., 0.;
} 

/******************************************************************************/
//
// Sets source strength
//
/******************************************************************************/
void Panel::setSourceStrength ( const double & sigin ) { _sigma = sigin; }

/******************************************************************************/
//
// Returns source strength
//
/******************************************************************************/
const double & Panel::sourceStrength () const { return _sigma; }

/******************************************************************************/
//
// Sets doublet strength
//
/******************************************************************************/
void Panel::setDoubletStrength ( const double & muin ) { _mu = muin; }

/******************************************************************************/
//
// Returns doublet strength
//
/******************************************************************************/
const double & Panel::doubletStrength () const { return _mu; }

/******************************************************************************/
//
// Returns panel area
//
/******************************************************************************/
const double & Panel::area () const { return _area; }

/******************************************************************************/
//
// Returns centroid
//
/******************************************************************************/
const Eigen::Vector3d & Panel::centroid () const { return _cen; }

/******************************************************************************/
//
// Returns normal vector
//
/******************************************************************************/
const Eigen::Vector3d & Panel::normal () const { return _norm; }

/******************************************************************************/
//
// Set / access flow velocity at centroid
//
/******************************************************************************/
void Panel::setVelocity ( const Eigen::Vector3d & vel ) { _vel = vel; }
const Eigen::Vector3d & Panel::velocity () const { return _vel; }
