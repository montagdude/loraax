#include <vector>
#include <string>
#include <Eigen/Core>
#include <cmath>
#include "util.h"
#include "algorithms.h"
#include "transformations.h"
#include "sectional_object.h"
#include "vertex.h"
#include "airfoil.h"
#include "section.h"

/******************************************************************************/
//
// Section class. Defines a wing section.
// 
/******************************************************************************/

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
Section::Section ()
{
  _xle = 0.;
  _zle = 0.;
  _chord = 0.;
  _twist = 0.;
  _roll = 0.;
  _nverts = 0;
  _verts.resize(0); 
  _uverts.resize(0); 
}

/******************************************************************************/
//
// Set or access position, orientation, and scale
//
/******************************************************************************/
void Section::setGeometry ( const double & xle, const double & y,
                            const double & zle, const double & chord,
                            const double & twist, const double & roll )
{
  _xle = xle;
  setY(y);
  _zle = zle;
  _chord = chord;
  _twist = twist;
  _roll = roll;
}

void Section::setRoll ( const double & roll ) { _roll = roll; }
const double & Section::xle () const { return _xle; }
const double & Section::zle () const { return _zle; }
const double & Section::chord () const { return _chord; }
const double & Section::twist () const { return _twist; }
const double & Section::roll () const { return _roll; }

/******************************************************************************/
//
// Access vertex
//
/******************************************************************************/
Vertex & Section::vert ( unsigned int idx )
{
  if ( idx >= _nverts )
  {
    conditional_stop(1, "Section::vert", "Index out of range.");
  }

  return _verts[idx];
} 

/******************************************************************************/
//
// Set vertices from spacing distribution
//   nchord: number of points along chord
//   lesprat: spacing at leading edge as fraction of uniform spacing
//   tesprat: spacing at leading edge as fraction of uniform spacing
//
/******************************************************************************/
void Section::setVertices ( unsigned int nchord, const double & lesprat,
                            const double & tesprat )
{
  Airfoil foil;
  double slen, sle, unisp, lesp, tesp;
  double a4top, a5top, a4bot, a5bot;
  std::vector<double> sv;
  std::vector<double> xf, zf;			// Vertices in foil coordinates
  unsigned int i;
  Eigen::Matrix3d rotation;

#ifdef DEBUG
  if (_foil.nBuffer() == 0)
  {
    conditional_stop(1, "Section::setVertices", "Airfoil not loaded yet.");
  }

  if (! _foil.splined())
  {
    conditional_stop(1, "Section::setVertices", "Must call splineFit first.");
  }
#endif

  // Save a copy of the section's airfoil and remove TE gap if necessary. By
  // working with a copy, the gap is only removed for 3D panel calculations, but
  // it is preserved for Xfoil BL calculations.

  foil = _foil;
  if (foil.teGap() > 1.E-14)
  {
    foil.modifyTEGap(0.0, 0.9);
    foil.smoothPaneling();
    foil.splineFit();
  }

  // Get spacings 

  slen = foil.sLen();
  sle = foil.sLE();
  unisp = slen / float(2*nchord-2);	// Uniform spacing
  lesp = unisp * lesprat;		// LE spacing
  tesp = unisp * tesprat;		// TE spacing

  // Optimize tanh spacing coefficients to minimize stretching

  opt_tanh_spacing(nchord, sle, tesp, lesp, a4top, a5top);
  opt_tanh_spacing(nchord, slen-sle, lesp, tesp, a4bot, a5bot);

  // Set spacing vector

  _nverts = 2*nchord - 1;
  _verts.resize(_nverts);
  _uverts.resize(_nverts);
  sv.resize(_nverts);
  sv[0] = 0.;
  for ( i = 1; i < nchord; i++ )
  {
    sv[i] = sv[i-1] +
      tanh_spacing(i-1, a4top, a5top, nchord, sle, tesp, lesp);
  }
  for ( i = 1; i < nchord; i++ )
  {
    sv[i+nchord-1] = sv[i+nchord-2] +
      tanh_spacing(i-1, a4bot, a5bot, nchord, slen-sle, lesp, tesp);
  }

  // Get vertices in foil coordinate system (unit chord, 0 <= x <= 1)

  xf.resize(_nverts);
  zf.resize(_nverts);
  for ( i = 0; i < _nverts; i++ )
  {
    foil.splineInterp(sv[i], xf[i], zf[i]);
  }

  // Transform to section coordinates. Also store non-rotated, non-translated
  // version for calculating sectional loads.

  rotation = inverse_euler_rotation(_roll, _twist, 0.0, "123");
  for ( i = 0; i < _nverts; i++ )
  {
    _verts[i].setCoordinates(xf[i], 0.0, zf[i]);
    _verts[i].translate(-0.25, 0., 0.);
    _verts[i].rotate(rotation);
    _verts[i].translate(0.25, 0., 0.);
    _verts[i].scale(_chord);
    _verts[i].translate(_xle, _y, _zle);

    _uverts[i].setCoordinates(xf[i], 0.0, zf[i]);
    _uverts[i].scale(_chord);
  } 
} 

/******************************************************************************/
//
// Access airfoil
//
/******************************************************************************/
Airfoil & Section::airfoil () { return _foil; }

/******************************************************************************/
//
// Computes pressure forces / span
//
/******************************************************************************/
void Section::computePressureForce ( const double & alpha, const double & uinf,
                                     const double & rhoinf )
{
  unsigned int i;
  double nx, nz, pave, qinf, lift, drag;
  Eigen::Vector3d force;
  Eigen::Matrix3d rotation;

  _fn = 0.;
  _fa = 0.;
  for ( i = 1; i < _nverts; i++ )
  {
    // Dimensional edge normal vector

    nx = _uverts[i].z() - _uverts[i-1].z();
    nz = _uverts[i-1].x() - _uverts[i].x();

    // Integrate pressure

    pave = 0.5*(_verts[i].data(5) + _verts[i-1].data(5));
    _fn -= pave*nz;
    _fa -= pave*nx;
  }

  // Rotate forces to inertial frame

  rotation = inverse_euler_rotation(_roll, _twist, 0.0, "123");
  force << _fa, 0., _fn;
  force = rotation*force;

  // Sectional lift and drag coefficients
 
  qinf = 0.5*rhoinf*std::pow(uinf,2.);
  lift = -force(0)*sin(alpha*M_PI/180.) + force(2)*cos(alpha*M_PI/180.);
  drag =  force(0)*cos(alpha*M_PI/180.) + force(2)*sin(alpha*M_PI/180.);
  _cl = lift/(qinf*_chord);
  _cd = drag/(qinf*_chord);
}

/******************************************************************************/
//
// Sectional lift and drag coefficients
//
/******************************************************************************/
const double & Section::liftCoefficient () const { return _cl; }
const double & Section::dragCoefficient () const { return _cd; }

/******************************************************************************/
//
// Computes Reynolds number based on local chord. Also sets is for airfoil.
//
/******************************************************************************/
void Section::computeReynoldsNumber ( const double & rhoinf,
                                     const double & uinf, const double & muinf )
{
  _re = rhoinf * uinf * _chord / muinf;
  _foil.setReynoldsNumber(_re);
} 

const double & Section::reynoldsNumber () const { return _re; }
