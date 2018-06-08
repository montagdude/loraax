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
// Interpolates BL data from airfoil points to section vertices and sets in
// vertex data. scale is optional and defaults to 1.
//
/******************************************************************************/
void Section::setVertexBLData ( const std::vector<double> & bldata,
                                unsigned int dataidx, double scale )
{
  unsigned int i, j1, j2;
  double interpval;

#ifdef DEBUG
  if ( int(bldata.size()) != _foil.nSmoothed() )
    conditional_stop(1, "Section::setVertexBLData", "Wrong size input vector.");
#endif

  for ( i = 0; i < _nverts; i++ )
  {
    j1 = _foilinterp[i].point1;
    j2 = _foilinterp[i].point2;
    interpval = bldata[j1]*_foilinterp[i].weight1
              + bldata[j2]*_foilinterp[i].weight2;
    _verts[i].setData(dataidx, interpval*scale);
  }
}

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
  _re = 0.;
  _converged = false;
  _foilinterp.resize(0);
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
  double sscale, svs;
  std::vector<double> sv;
  std::vector<double> xf, zf;			// Vertices in foil coordinates
  unsigned int i, j, nsmoothed;
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

  // Finds interpolation points on airfoil for section vertices

  _foilinterp.resize(_nverts);
  nsmoothed = _foil.nSmoothed();

  sscale = _foil.sSmoothed(nsmoothed-1) / (sv[_nverts-1] + 1.E-12);
    // Scale section to airfoil length, leaving a little room for roundoff error

  for ( i = 0; i < _nverts; i++ )
  {
    _foilinterp[i].weight1 = 0.;
    _foilinterp[i].weight2 = 0.;
    svs = sv[i] * sscale;
    for ( j = 0; j < nsmoothed-1; j++ )
    {
      if ( (svs >= _foil.sSmoothed(j)) && (svs <= _foil.sSmoothed(j+1)) )
      {
        _foilinterp[i].point1 = j;
        _foilinterp[i].point2 = j+1;
        _foilinterp[i].weight2 = (svs - _foil.sSmoothed(j))
                               / (_foil.sSmoothed(j+1) - _foil.sSmoothed(j));
        _foilinterp[i].weight1 = 1. - _foilinterp[i].weight2;
        break;
      }
    }
#ifdef DEBUG
    if ( (_foilinterp[i].weight1 == 0.) && (_foilinterp[i].weight2 == 0.) )
      print_warning("Section::setVertices",
                    "Could not find airfoil->section interpolants." );
#endif
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

/******************************************************************************/
//
// Sets Mach number for airfoil
//
/******************************************************************************/
void Section::setMachNumber ( const double & minf )
{
  _foil.setMachNumber(minf); 
}

/******************************************************************************/
//
// Boundary layer calculations with Xfoil
//
/******************************************************************************/
void Section::computeBL ( const Eigen::Vector3d & uinfvec,
                          const double & rhoinf, const double & alpha )
{
  Eigen::Vector3d uinfvec_p;
  double qinf, uinf, uinfp, cl2d;
  Eigen::Matrix3d inertial2section;
  std::vector<double> bldata;
  int stat;

  qinf = 0.5*rhoinf*uinfvec.squaredNorm();
  uinf = uinfvec.norm();

  // Sectional lift must be computed as an input to Xfoil. Sectional forces and
  // moments will be recomputed after running Xfoil for the purpose of writing
  // data to the sectional output files.

  computeForceMoment(alpha, uinf, rhoinf, true);

  /** To get 2D Cl:
      1. Transform uinfvec to section frame -> uinfvec_p
      2. 2D angle of attack is atan(uinfvec_p[2]/uinfvec_p[0])
      3. 2D lift is -_fa*sin(alpha2d) + _fn*cos(alpha2d)
  **/
  //FIXME: this needs to be checked

  inertial2section = euler_rotation(_roll, _twist, 0.0, "123");
  uinfvec_p = inertial2section.transpose() * uinfvec;
  uinfp = uinfvec_p.norm();
  cl2d = -_fa*uinfvec_p[2]/uinfp + _fn*uinfvec_p[0]/uinfp;
  cl2d /= qinf*_chord;

  // Run xfoil at 2D Cl

  if (_foil.runXfoil(cl2d) != 0)
    _converged = false;
  else
    _converged = true;

  // Interpolate BL quantities to vertices. These will be overwritten for
  // unconverged sections if interpolation/extrapolation is possible. Perform
  // scaling as needed.

  bldata = _foil.blData("cf", stat);
  setVertexBLData(bldata, 7);
  bldata = _foil.blData("deltastar", stat);
  setVertexBLData(bldata, 8, _chord);
  bldata = _foil.blData("ampl", stat);
  setVertexBLData(bldata, 9);
}

bool Section::blConverged () const { return _converged; }

/******************************************************************************/
//
// Computes force and moment / unit span
//
/******************************************************************************/
void Section::computeForceMoment ( const double & alpha, const double & uinf,
                                   const double & rhoinf, bool viscous )
{
  unsigned int i;
  double nx, nz, tx, tz, cenx, cenz, pave, cfave, qinf;
  double fap, fav, fnp, fnv, mp, mv, liftp, liftv, dragp, dragv;
  double dfap, dfav, dfnp, dfnv;
  Eigen::Vector3d forcep, forcev, momentp, momentv;
  Eigen::Matrix3d section2inertial;

  fap = 0.;
  fav = 0.;
  fnp = 0.;
  fnv = 0.;
  mp = 0.;
  mv = 0.;
  qinf = 0.5*rhoinf*std::pow(uinf,2.);
  for ( i = 1; i < _nverts; i++ )
  {
    // Dimensional edge normal vector and edge center

    nx = _uverts[i].z() - _uverts[i-1].z();
    nz = _uverts[i-1].x() - _uverts[i].x();
	cenx = 0.5*(_uverts[i].x() + _uverts[i-1].x());
	cenz = 0.5*(_uverts[i].z() + _uverts[i-1].z());

    // Integrate pressure

    pave = 0.5*(_verts[i].data(5) + _verts[i-1].data(5));
    dfnp = -pave*nz;
    dfap = -pave*nx;
	fnp += dfnp;
	fap += dfap;
	mp += dfnp * (0.25*_chord - cenx) + dfap*cenz;

	// Viscous contribution

	if (viscous)
	{
	  tx = _uverts[i].x() - _uverts[i-1].x();
	  tz = _uverts[i].z() - _uverts[i-1].z();
	  if (tx < 0)
	  {
		tx *= -1.;
		tz *= -1.;
	  }
      cfave = 0.5*(_verts[i].data(7) + _verts[i-1].data(7));
	  dfnv = cfave*qinf*tz;
	  dfav = cfave*qinf*tx;
	  fnv += dfnv;
	  fav += dfav;
	  mv += dfnv * (0.25*_chord - cenx) + dfav*cenz;
	}
  }
  _fn = fnp + fnv;
  _fa = fap + fav;

  // Rotate forces and moments to inertial frame

  section2inertial = inverse_euler_rotation(_roll, _twist, 0.0, "123");
  forcep << fap, 0., fnp;
  forcep = section2inertial*forcep;
  forcev << fav, 0., fnv;
  forcev = section2inertial*forcev;
  momentp << 0., mp, 0.;
  momentp = section2inertial*momentp;
  momentv << 0., mv, 0.;
  momentv = section2inertial*momentv;

  // Sectional lift, drag, and moment coefficients
 
  liftp = -forcep(0)*sin(alpha*M_PI/180.) + forcep(2)*cos(alpha*M_PI/180.);
  liftv = -forcev(0)*sin(alpha*M_PI/180.) + forcev(2)*cos(alpha*M_PI/180.);
  dragp =  forcep(0)*cos(alpha*M_PI/180.) + forcep(2)*sin(alpha*M_PI/180.);
  dragv =  forcev(0)*cos(alpha*M_PI/180.) + forcev(2)*sin(alpha*M_PI/180.);

  _clp = liftp/(qinf*_chord);
  _clv = liftv/(qinf*_chord);
  _cdp = dragp/(qinf*_chord);
  _cdv = dragv/(qinf*_chord);
  _cmp = momentp(1)/(qinf*_chord*_chord);
  _cmv = momentv(1)/(qinf*_chord*_chord);
}

/******************************************************************************/
//
// Sectional force and moment coefficients
//
/******************************************************************************/
double Section::liftCoefficient () const { return _clp + _clv; }
const double & Section::pressureLiftCoefficient () const { return _clp; }
const double & Section::viscousLiftCoefficient () const { return _clv; }

double Section::dragCoefficient () const { return _cdp + _cdv; }
const double & Section::pressureDragCoefficient () const { return _cdp; }
const double & Section::viscousDragCoefficient () const { return _cdv; }

double Section::pitchingMomentCoefficient () const { return _cmp + _cmv; }
const double & Section::pressurePitchingMomentCoefficient () const { return _cmp; }
const double & Section::viscousPitchingMomentCoefficient () const { return _cmv; }
