#include <fstream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "sectional_object.h"
extern "C"
{
  #include <xfoil_interface.h>
}
#include "util.h"
#include "airfoil.h"

/******************************************************************************/
//
// Airfoil class. Stores and manipulates airfoil coordinates and solves BL
// equations with Xfoil.
//
/******************************************************************************/

/******************************************************************************/
//
// Helper function to avoid reusing code between copy constructor and copy
// assignment
//
/******************************************************************************/
void Airfoil::copyData ( const Airfoil & foil )
{
  xfoil_copy(&foil._xdg, &_xdg);
  _nb = foil._nb;
  _n = foil._n;
  _s = foil._s;
  _xs = foil._xs;
  _zs = foil._zs;
  _sle = foil._sle;
  _unit_transform = foil._unit_transform;
  _cp = foil._cp;
  _cf = foil._cf;
}

/******************************************************************************/
//
// Constructor, copy constructor, copy assignment, destructor
//
/******************************************************************************/
Airfoil::Airfoil ()
{
  xfoil_init(&_xdg);
  _nb = 0;
  _n = 0;
  _s.resize(0);
  _xs.resize(0);
  _zs.resize(0);
  _y = 0.;
  _unit_transform = false;
  _cp.resize(0);
  _cf.resize(0);
}

Airfoil::Airfoil ( const Airfoil & foil )
{
  xfoil_init(&_xdg);
  copyData(foil);
}

Airfoil & Airfoil::operator = ( const Airfoil & foil )
{
  copyData(foil);
  return *this;
}

Airfoil::~Airfoil()
{
  xfoil_cleanup(&_xdg);
  _nb = 0;
  _n = 0;
  _s.resize(0);
  _xs.resize(0);
  _zs.resize(0);
  _y = 0.;
  _unit_transform = false;
  _cp.resize(0);
  _cf.resize(0);
}

/******************************************************************************/
//
// Reads airfoil coordinates from file. Returns 0 on success, 1 on failure to
// open file, 2 on format error.
//
/******************************************************************************/
int Airfoil::readCoordinates ( const std::string & fname )
{
  std::ifstream f;
  std::string line;
  std::vector<double> xb, zb;
  std::vector<std::string> xz;
  double x, z;
  int i;
  int retval = 0;

  f.open(fname.c_str());
  if (not f.is_open()) { return 1; }

  // Start reading coordinates or skip if header line

  _nb = 0;
  xb.resize(0);
  zb.resize(0);
  std::getline(f, line);
  xz = split_string(bracket_name(line));
  if ( (xz.size() == 2) && (string2double(xz[0], x) == 0) &&
       (string2double(xz[1], z) == 0) )
  {
    _nb += 1;
    xb.push_back(x);
    zb.push_back(z);
  }

  // Read rest of coordinates

  while (1)
  {
    std::getline(f, line);
    if (f.eof()) { break; }

    xz = split_string(bracket_name(line));
    if ( (xz.size() != 2) || (string2double(xz[0], x) != 0) ||
         (string2double(xz[1], z) != 0) )	// Read error
    {
      retval = 2;
      break;
    }
    else
    {
      _nb += 1;
      xb.push_back(x);
      zb.push_back(z);
    }
  }
  f.close();

  // Set coordinates in Xfoil

  if (retval == 0)
  {
    double xba[_nb], zba[_nb];
    for ( i = 0; i < _nb; i++ )
    {
      xba[i] = xb[i];
      zba[i] = zb[i];
    }
    xfoil_set_buffer_airfoil(&_xdg, xba, zba, &_nb);
  }

  _unit_transform = false;

  return retval;
}

/******************************************************************************/
//
// Sets coordinates based on 4-digit NACA designation
//
/******************************************************************************/
int Airfoil::naca4Coordinates ( const std::string & des,
                                const int & npointside )
{
  double x[2*npointside], z[2*npointside];
  int i;

  if (des.size() != 4)
    return 1;

  naca_4_digit(des.c_str(), &npointside, x, z, &_nb);
  double xba[_nb], zba[_nb];
  for ( i = 0; i < _nb; i++ )
  {
    xba[i] = x[i];
    zba[i] = z[i];
  }
  xfoil_set_buffer_airfoil(&_xdg, xba, zba, &_nb);

  _unit_transform = false;

  return 0;
}

/******************************************************************************/
//
// Sets coordinates based on 5-digit NACA designation
//
/******************************************************************************/
int Airfoil::naca5Coordinates ( const std::string & des,
                                const int & npointside )
{
  double x[2*npointside], z[2*npointside];
  int i, stat;

  if (des.size() != 5)
    return 1;

  naca_5_digit(des.c_str(), &npointside, x, z, &_nb, &stat);
  if (stat != 0)
    return 2;

  double xba[_nb], zba[_nb];
  for ( i = 0; i < _nb; i++ )
  {
    xba[i] = x[i];
    zba[i] = z[i];
  }
  xfoil_set_buffer_airfoil(&_xdg, xba, zba, &_nb);

  _unit_transform = false;

  return 0;
}

/******************************************************************************/
//
// Sets airfoil coordinates from given arrays
//
/******************************************************************************/
int Airfoil::setCoordinates ( const std::vector<double> & x,
                              const std::vector<double> & z )
{
  int i;

  if (x.size() != z.size())
    return 1;

  _nb = x.size();
  double xba[_nb], zba[_nb];
  for ( i = 0; i < _nb; i++ )
  {
    xba[i] = x[i];
    zba[i] = z[i];
  }
  xfoil_set_buffer_airfoil(&_xdg, xba, zba, &_nb);

  _unit_transform = false;

  return 0;
}

/******************************************************************************/
//
// Sets airfoil buffer coordinates by interpolating smoothed coordinates from
// two other airfoils. It is assumed that both airfoils have the same number
// of vertices in smoothed coordinates arrays, whereas buffer coordinates may
// differ.
//
// interpfrac: 0 to 1. If 0, output = foil1 coordinates; if 1, output = foil2
//             coordinates.
//
/******************************************************************************/
int Airfoil::interpCoordinates ( const Airfoil & foil1, const Airfoil & foil2,
                                 const double interpfrac )
{
  int i;
  std::vector<double> x1, z1, x2, z2;

  if ( (foil1.nSmoothed() == 0) || (foil2.nSmoothed() == 0) )
  {
#ifdef DEBUG
    print_warning("Airfoil::interpCoordinates",
                  "Smoothed coordinates not present for foil1 or foil2.");
#endif
    return 1;
  }

  if (foil1.nSmoothed() != foil2.nSmoothed() )
  {
#ifdef DEBUG
    print_warning("Airfoil::interpCoordinates",
                  "foil1 and foil2 must have same number of smoothed panels");
#endif
    return 2;
  }

  _nb = foil1.nSmoothed();
  double xba[_nb], zba[_nb];
  foil1.smoothedCoordinates(x1, z1);
  foil2.smoothedCoordinates(x2, z2);
  for ( i = 0; i < _nb; i++ )
  {
    xba[i] = (1.-interpfrac)*x1[i] + interpfrac*x2[i];
    zba[i] = (1.-interpfrac)*z1[i] + interpfrac*z2[i];
  }
  xfoil_set_buffer_airfoil(&_xdg, xba, zba, &_nb);

  _unit_transform = false;

  return 0;
}

/******************************************************************************/
//
// Puts buffer coordinates in counterclockwise order. Does nothing if already
// in ccw order.
//
/******************************************************************************/
void Airfoil::ccwOrderCoordinates ()
{
  double cross;
  int i;
  double dx1, dx2, dz1, dz2;
  double xrev[_nb], zrev[_nb];

#ifdef DEBUG
  if (_nb == 0)
  {
    conditional_stop(1, "Airfoil::ccwOrderCoordinates",
                     "Airfoil not yet loaded.");
  }
#endif

  // Cross product e2 x e1 integral will be in -y direction for ccw ordering

  cross = 0.;
  for ( i = 1; i < _nb-1; i++ )
  {
    dx1 = _xdg.xfd.XB[i-1] - _xdg.xfd.XB[i];
    dx2 = _xdg.xfd.XB[i+1] - _xdg.xfd.XB[i];
    dz1 = _xdg.xfd.YB[i-1] - _xdg.xfd.YB[i];
    dz2 = _xdg.xfd.YB[i+1] - _xdg.xfd.YB[i];
    cross += (dz2*dx1 - dx2*dz1);
  }
  if (cross < 0.) { return; }

  // Reverse ordering for clockwise airfoil

  for ( i = 0; i < _nb; i++ )
  {
    xrev[i] = _xdg.xfd.XB[_nb-i-1];
    zrev[i] = _xdg.xfd.YB[_nb-i-1];
  }
  xfoil_set_buffer_airfoil(&_xdg, xrev, zrev, &_nb);
}

/******************************************************************************/
//
// Fits a spline to buffer coordinates
//
/******************************************************************************/
int Airfoil::splineFit ()
{
  double s[_nb], xs[_nb], zs[_nb];
  int i;

  if (_nb == 0)
  {
#ifdef DEBUG
    print_warning("Airfoil::splineFit", "Airfoil not yet loaded.");
#endif
    return 1;
  }

  _s.resize(_nb);
  _xs.resize(_nb);
  _zs.resize(_nb);
  xfoil_spline_coordinates(_xdg.xfd.XB, _xdg.xfd.YB, &_nb, s, xs, zs);

  for ( i = 0; i < _nb; i++ )
  {
    _s[i] = s[i];
    _xs[i] = xs[i];
    _zs[i] = zs[i];
  }

  return 0;
}

/******************************************************************************/
//
// Retrieve spline information
//
/******************************************************************************/
const double & Airfoil::sLen () const { return _s.back(); }
const double & Airfoil::sLE () const { return _sle; }
bool Airfoil::splined () const
{
  if (_s.size() > 0)
    return true;
  else
    return false;
}

/******************************************************************************/
//
// Compute coordinate at spline parameter
//
/******************************************************************************/
int Airfoil::splineInterp ( const double & sc, double & xc, double & zc ) const
{
  double s[_nb], xs[_nb], zs[_nb];
  int i;

  if (! splined())
  {
#ifdef DEBUG
    print_warning("Airfoil::splineInterp", "Must call splineFit first.");
#endif
    return 1;
  }

  for ( i = 0; i < _nb; i++ )
  {
    s[i] = _s[i];
    xs[i] = _xs[i];
    zs[i] = _zs[i];
  }
  xfoil_eval_spline(_xdg.xfd.XB, _xdg.xfd.YB, s, xs, zs, &_nb, &sc, &xc, &zc);

  return 0;
}

/******************************************************************************/
//
// Scales to unit chord and moves to origin.
//
/******************************************************************************/
int Airfoil::unitTransform ()
{
  double xle, zle, xte, chord;
  double s[_nb], xs[_nb], zs[_nb];
  int i;

  if (_nb == 0)
  {
#ifdef DEBUG
    print_warning("Airfoil::unitTransform", "Airfoil not yet loaded.");
#endif
    return 1;
  }

  if (_s.size() == 0)
  {
#ifdef DEBUG
    print_warning("Airfoil::unitTransform", "Must call splineFit first.");
#endif
    return 2;
  }

  for ( i = 0; i < _nb; i++ )
  {
    s[i] = _s[i];
    xs[i] = _xs[i];
    zs[i] = _zs[i];
  }
  xfoil_lefind(_xdg.xfd.XB, _xdg.xfd.YB, s, xs, zs, &_nb, &_sle, &xle, &zle);

  xte = -1E+06;
  for ( i = 0; i < _nb; i++ )
  {
    if (_xdg.xfd.XB[i] > xte) { xte = _xdg.xfd.XB[i]; }
  }
  chord = xte - xle;
  if (chord <= 0.)
  {
#ifdef DEBUG
    print_warning("Airfoil::unitTransform", "Airfoil has chord <= 0.");
#endif
    return 3;
  }

  // Move leading edge to origin and scale

  for ( i = 0; i < _nb; i++ )
  {
    _xdg.xfd.XB[i] -= xle;
    _xdg.xfd.YB[i] -= zle;
    _xdg.xfd.XB[i] /= chord;
    _xdg.xfd.YB[i] /= chord;
  }

  // Refit spline to transformed coordinates

  splineFit();
  for ( i = 0; i < _nb; i++ )
  {
    s[i] = _s[i];
    xs[i] = _xs[i];
    zs[i] = _zs[i];
  }
  xfoil_lefind(_xdg.xfd.XB, _xdg.xfd.YB, s, xs, zs, &_nb, &_sle, &xle, &zle);

  _unit_transform = true;

  return 0;
}

/******************************************************************************/
//
// Sets Xfoil paneling and run options
//
/******************************************************************************/
void Airfoil::setXfoilOptions ( const xfoil_options_type & xfoil_opts,
                                const xfoil_geom_options_type & geom_opts )
{
  xfoil_defaults(&_xdg, &xfoil_opts);
  xfoil_set_paneling(&_xdg, &geom_opts);
  _n = geom_opts.npan;
}

/******************************************************************************/
//
// Generates smoothed coordinates. Returns 1 and does nothing if buffer airfoil
// has not been loaded. Returns 2 and does nothing if unitTransform has not been
// called. Returns 3 if smoothing failed.
//
/******************************************************************************/
int Airfoil::smoothPaneling ()
{
  int stat;

  if (_nb == 0)
  {
#ifdef DEBUG
    print_warning("Airfoil::smoothPaneling", "Airfoil not yet loaded.");
#endif
    return 1;
  }

  if (! _unit_transform)
  {
#ifdef DEBUG
    print_warning("Airfoil::smoothPaneling", "Must call unitTransform first.");
#endif
    return 2;
  }

  xfoil_smooth_paneling(&_xdg, &stat);
  if (stat != 0)
  {
    conditional_stop(1, "Airfoil::smoothPaneling", "Airfoil smoothing failed.");
    return 3;
  }

  return 0;
}

/******************************************************************************/
//
// Returns TE gap
//
/******************************************************************************/
double Airfoil::teGap () const
{
  double dx, dz, gap;

  dx = _xdg.xfd.XB[0] - _xdg.xfd.XB[_nb-1];
  dz = _xdg.xfd.YB[0] - _xdg.xfd.YB[_nb-1];
  gap = std::sqrt(std::pow(dx,2.) + std::pow(dz,2.));

  return gap;
}

/******************************************************************************/
//
// Modifies trailing edge gap. Returns 1 and does nothing if buffer airfoil has
// not been loaded. Returns 2 and does nothing if unitTransform has not been
// called. Returns 3 if TE gap modification failed.
//
// Note: this method only modifies the buffer coordinates. You may want to call
// splineInterp, smoothPaneling, etc. again afterwards.
//
/******************************************************************************/
int Airfoil::modifyTEGap ( const double & newgap, const double & blendloc )
{
  int stat;

  if (_nb == 0)
  {
#ifdef DEBUG
    print_warning("Airfoil::modifyTEGap", "Airfoil not yet loaded.");
#endif
    return 1;
  }

  if (! _unit_transform)
  {
#ifdef DEBUG
    print_warning("Airfoil::modifyTEGap", "Must call unitTransform first.");
#endif
    return 2;
  }

  xfoil_modify_tegap(&_xdg, &newgap, &blendloc, &_nb, &stat);
  if (stat != 0)
  {
    conditional_stop(1, "Airfoil::modifyTEGap", "TE gap modification failed.");
    return 3;
  }

  return 0;
}

/******************************************************************************/
//
// Get airfoil data
//
/******************************************************************************/
int Airfoil::nBuffer () const { return _nb; }
int Airfoil::nSmoothed () const { return _n; }
void Airfoil::bufferCoordinates ( std::vector<double> & xb,
                                  std::vector<double> & zb ) const
{
  int i;

#ifdef DEBUG
  if (_nb == 0)
  {
    print_warning("Airfoil::bufferCoordinates",
                  "Airfoil not yet loaded.");
  }
#endif

  xb.resize(_nb);
  zb.resize(_nb);
  for ( i = 0; i < _nb; i++ )
  {
    xb[i] = _xdg.xfd.XB[i];
    zb[i] = _xdg.xfd.YB[i];
  }
}

void Airfoil::smoothedCoordinates ( std::vector<double> & x,
                                    std::vector<double> & z ) const
{
  int i;

#ifdef DEBUG
  if (_n == 0)
  {
    print_warning("Airfoil::smoothedCoordinates",
                  "Smoothed coordinates not yet available.");
  }
#endif

  x.resize(_n);
  z.resize(_n);
  for ( i = 0; i < _n; i++ )
  {
    x[i] = _xdg.xfd.X[i];
    z[i] = _xdg.xfd.Y[i];
  }
}

/******************************************************************************/
//
// Set Reynolds number in Xfoil
//
/******************************************************************************/
void Airfoil::setReynoldsNumber ( const double & re )
{
  xfoil_set_reynolds_number(&_xdg, &re);
}

/******************************************************************************/
//
// Set Mach number in Xfoil
//
/******************************************************************************/
void Airfoil::setMachNumber ( const double & mach )
{
  xfoil_set_mach_number(&_xdg, &mach);
}

/******************************************************************************/
//
// Run Xfoil at specified lift coefficient
// stat: 0 for success, 1 if convergence failed
//
/******************************************************************************/
int Airfoil::runXfoil ( const double & clspec )
{
  double alpha, cl, cd, cm;
  bool converged;
  int stat;
  
  xfoil_speccl(&_xdg, &clspec, &alpha, &cl, &cd, &cm, &converged, &stat);
  if (converged)
    return 0;
  else
    return 1;
}
