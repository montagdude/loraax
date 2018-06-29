#include <vector>
#include <string>
#include <cmath>
#include <Eigen/Dense>
#include "util.h"
#include "algorithms.h"
#include "vertex.h"
#include "element.h"
#include "panel.h"

/******************************************************************************/
//
// Panel class. Computes geometric quantities, source/doublet influence
// coefficients at a point, etc.
// 
/******************************************************************************/

const double Panel::_farfield_distance_factor = 8.;

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
	_norm << 0., 0., 0.;
	_tan << 0., 0., 0.;
	_xtrans.resize(0);
	_vel << 0., 0., 0.;
	_p = 0.;
	_cp = 0.;
	_cf = 0.;
	_mdefect = 0.;
	_dmdefect = 0.;
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
// Computes and performs LU factorization of grid transformation. Also computes
// surface tangential direction for viscous forces.
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
// Sets source strength to provided value
//
/******************************************************************************/
void Panel::setSourceStrength ( const double & sigma ) { _sigma = sigma; }

/******************************************************************************/
//
// Computes and sets source strength
//
/******************************************************************************/
void Panel::computeSourceStrength ( const Eigen::Vector3d & uinfvec,
                                    bool viscous )
{
	double dsp, dsm, mdefp, mdefm;

	// Compute mass defect derivative for viscous cases

	if (viscous)
	{
		if (_front == NULL)
		{
			if (_back == NULL)
			{
				conditional_stop(1, "Panel::computeSourceStrength",
				                "Must have at least a front or back neighbor.");
			}
			dsm = -(_cen - _back->centroid()).norm();
			mdefm = _back->massDefect();
			_dmdefect = -(_mdefect - mdefm) / dsm;
		}
		else if (_back == NULL)
		{
			dsp = (_cen - _front->centroid()).norm();
			mdefp = _front->massDefect();
			_dmdefect = (mdefp - _mdefect) / dsp;
		}
		else
		{
			dsp = (_cen - _front->centroid()).norm();
			mdefp = _front->massDefect();
			dsm = -(_cen - _back->centroid()).norm();
			mdefm = _back->massDefect();
			_dmdefect = centered_difference(mdefp, _mdefect, mdefm, dsp, dsm);
		}

		// Set correct sign based on edge velocity

		if (_vel.transpose() * _tan < 0.)
			_dmdefect *= -1.;
	}

	// Source strength

	_sigma = -uinfvec.transpose()*_norm + _dmdefect;
}

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
// Returns surface tangent vector
//
/******************************************************************************/
const Eigen::Vector3d & Panel::tan () const { return _tan; }

/******************************************************************************/
//
// Set usrface tangent vector
//
/******************************************************************************/
void Panel::setTangent ( const Eigen::Vector3d & tan ) { _tan = tan; }

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
// Mass defect, uedge*deltastar, and surface derivative
//
/******************************************************************************/
const double & Panel::massDefect () const { return _mdefect; }
const double & Panel::massDefectDerivative () const { return _dmdefect; }

/******************************************************************************/
//
// Averages viscous quantities from vertices
//
/******************************************************************************/
void Panel::averageFromVertices ()
{
	unsigned int i, nverts;
	double dx, dy, dz, dist, weightsum;

	nverts = _verts.size();
	_cf = 0.; 
	_mdefect = 0.; 
	weightsum = 0.;
	for ( i = 0; i < nverts; i++ )
	{
		dx = _verts[i]->x() - _cen(0);
		dy = _verts[i]->y() - _cen(1);
		dz = _verts[i]->z() - _cen(2);
		dist = std::sqrt(std::pow(dx,2.) + std::pow(dy,2.) +
		                 std::pow(dz,2.));
		_cf += _verts[i]->data(7)/dist;
		_mdefect += (_verts[i]->data(8)*_verts[i]->data(10))/dist;
		weightsum += 1./dist;
	}
	_cf /= weightsum;
	_mdefect /= weightsum;
}

/******************************************************************************/
//
// Compute force and moment contributions
//
/******************************************************************************/
void Panel::computeForceMoment ( const double & uinf, const double & rhoinf,
                                 const double & pinf,
                                 const Eigen::Vector3d & moment_center,
                                 bool viscous, Eigen::Vector3d & fp,
                                 Eigen::Vector3d & fv, Eigen::Vector3d & mp,
                                 Eigen::Vector3d & mv )
{
	double q, tau, sign;

	// Pressure force and moment

	q = 0.5*rhoinf*std::pow(uinf,2.);
	fp = -(_p-pinf)*_norm*_area;
	mp = (_cen - moment_center).cross(fp);

	if (viscous)
	{
		// Average viscous quantities from vertices

		averageFromVertices();

		// Determine sign on _tan vector corresponding with flow direction

		if (_vel.transpose()*_tan < 0.)
			sign = -1.;
		else
			sign = 1.;

		// Compute viscous force and moment

		tau = _cf * q;
		fv = sign*tau*_tan*_area;
		mv = (_cen - moment_center).cross(fv);
	}
	else
	{
		fv << 0., 0., 0.;
		mv << 0., 0., 0.;
	}
}
