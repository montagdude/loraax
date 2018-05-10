#include <vector>
#include <string>
#include <cmath>
#include <Eigen/Dense>
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
  _cp = 0.;
  _p = 0.;
  _colloc << 0., 0., 0.;
  _colloc_is_centroid = true;
  _right = NULL;
  _left = NULL;
  _front = NULL;
  _back = NULL;
} 

/******************************************************************************/
//
// Set/access neighbor panels
//
/******************************************************************************/
void Panel::setRightNeighbor ( Panel * right )
{
  if (_right != NULL)
  {
    print_warning("Panel::setRightNeighbor",
                  "Resetting right neighbor panel.");
  }
  _right = right;
}

void Panel::setLeftNeighbor ( Panel * left )
{
  if (_left != NULL)
  {
    print_warning("Panel::setLeftNeighbor",
                  "Resetting left neighbor panel.");
  }
  _left = left;
}

void Panel::setFrontNeighbor ( Panel * front )
{
  if (_front != NULL)
  {
    print_warning("Panel::setFrontNeighbor",
                  "Resetting front neighbor panel.");
  }
  _front = front;
}

void Panel::setBackNeighbor ( Panel * back )
{
  if (_back != NULL)
  {
    print_warning("Panel::setBackNeighbor",
                  "Resetting back neighbor panel.");
  }
  _back = back;
}

Panel * Panel::rightNeighbor () { return _right; }
Panel * Panel::leftNeighbor () { return _left; }
Panel * Panel::frontNeighbor () { return _front; }
Panel * Panel::backNeighbor () { return _back; }

/******************************************************************************/
//
// Computes and performs LU factorization of grid transformation
//
/******************************************************************************/
int Panel::computeGridTransformation ()
{
  double dxdxi, dydxi, dzdxi, dxdeta, dydeta, dzdeta, dxdchi, dydchi, dzdchi;

  if (_front == NULL)
  {
    if (_back == NULL)
    {
      conditional_stop(1, "Panel::computeGridTransformation",
                       "Must have at least a front or back neighbor."); 
      return 1;
    }
    dxdeta = _cen(0) - _back->centroid()(0);
    dydeta = _cen(1) - _back->centroid()(1);
    dzdeta = _cen(2) - _back->centroid()(2);
  }
  else if (_back == NULL)
  {
    dxdeta = _front->centroid()(0) - _cen(0);
    dydeta = _front->centroid()(1) - _cen(1);
    dzdeta = _front->centroid()(2) - _cen(2);
  }
  else
  {
    dxdeta = 0.5*(_front->centroid()(0) - _back->centroid()(0));
    dydeta = 0.5*(_front->centroid()(1) - _back->centroid()(1));
    dzdeta = 0.5*(_front->centroid()(2) - _back->centroid()(2));
  }

  if (_right == NULL)
  {
    if (_left == NULL)
    {
      conditional_stop(1, "Panel::computeGridTransformation",
                       "Must have at least a right or left neighbor."); 
      return 1;
    }
    dxdxi = _cen(0) - _left->centroid()(0);
    dydxi = _cen(1) - _left->centroid()(1);
    dzdxi = _cen(2) - _left->centroid()(2);
  }
  else if (_left == NULL)
  {
    dxdxi = _right->centroid()(0) - _cen(0);
    dydxi = _right->centroid()(1) - _cen(1);
    dzdxi = _right->centroid()(2) - _cen(2);
  }
  else
  {
    dxdxi = 0.5*(_right->centroid()(0) - _left->centroid()(0));
    dydxi = 0.5*(_right->centroid()(1) - _left->centroid()(1));
    dzdxi = 0.5*(_right->centroid()(2) - _left->centroid()(2));
  }

  dxdchi = _norm(0);
  dydchi = _norm(1);
  dzdchi = _norm(2);

  // Grid jacobian matrix and factorization

  _jac << dxdxi,  dydxi,  dzdxi,
          dxdeta, dydeta, dzdeta,
          dxdchi, dydchi, dzdchi;
  _lu.compute(_jac);

  return 0;
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
// Set/access collocation point
//
/******************************************************************************/
void Panel::setCollocationPoint ( const Eigen::Vector3d & colloc )
{
  _colloc = colloc;
  _colloc_is_centroid = false;
}

const Eigen::Vector3d & Panel::collocationPoint () const
{
  if (_colloc_is_centroid)
    return _cen;
  else
    return _colloc;
}

bool Panel::collocationPointIsCentroid () const { return _colloc_is_centroid; }

/******************************************************************************/
//
// Compute / access flow velocity at centroid
//
/******************************************************************************/
void Panel::computeVelocity ( const Eigen::Vector3d & uinfvec )
{
  double dmudxi, dmudeta, dmudchi;
  Eigen::Vector3d gradmu_grid, gradmu;

  // Compute surface gradient of doublet strength

  if (_front == NULL)
    dmudeta = _mu - _back->doubletStrength();
  else if (_back == NULL)
    dmudeta = _front->doubletStrength() - _mu;
  else
    dmudeta = 0.5*(_front->doubletStrength() - _back->doubletStrength());

  if (_right == NULL)
    dmudxi = _mu - _left->doubletStrength();
  else if (_left == NULL)
    dmudxi = _right->doubletStrength() - _mu;
  else
    dmudxi = 0.5*(_right->doubletStrength() - _left->doubletStrength());

  dmudchi = 0.;

  gradmu_grid << dmudxi, dmudeta, dmudchi;
  gradmu = _lu.solve(gradmu_grid);

  // Total velocity

  _vel = gradmu + _sigma*_norm + uinfvec; 
}
const Eigen::Vector3d & Panel::velocity () const { return _vel; }

/******************************************************************************/
//
// Compute / access pressure and pressure coefficient. Also performs
// compressibility corrections on surface pressure and velocity.
//
/******************************************************************************/
int Panel::computePressure ( const double & uinf, const double & rhoinf,
                             const double & pinf, const double & minf,
                             bool compressible )
{
  double uinf2, vel2, cpinc, minf2, beta, lambda, den;

  uinf2 = std::pow(uinf, 2.);
  vel2 = _vel.squaredNorm();
  cpinc = 1. - vel2 / uinf2;

  // Karman-Tsien compressibility correction

  if (compressible)
  {
    minf2 = std::pow(minf, 2.);
    beta = std::sqrt(1. - minf2);
    lambda = minf2 / std::pow(1 + beta, 2.);
    den = (beta+lambda*(1+beta)*0.5*cpinc); 
    if (den <= 0.0)
    {
      conditional_stop(1, "Panel::computePressure",
                       "Locally supersonic flow detected. " +
                       std::string("Compressibility correction is invalid."));
      return 1;
    }
    _cp = cpinc / den;
    _vel *= (1. - lambda) / (1. - lambda*vel2/uinf2);
  }
  else
    _cp = cpinc;

  _p = 0.5*_cp*rhoinf*uinf2 + pinf;

  return 0;
}

const double & Panel::pressure () const { return _p; }
const double & Panel::pressureCoefficient () const { return _cp; }

/******************************************************************************/
//
// Compute force and moment contributions
//
/******************************************************************************/
void Panel::computeForceMoment ( const double & pinf, Eigen::Vector3d & force,
                                 Eigen::Vector3d & moment,
                                 const Eigen::Vector3d & moment_center ) const
{
//FIXME: add viscous terms
  force = -(_p-pinf)*_norm*_area;
  moment = (_cen - moment_center).cross(force);
} 
