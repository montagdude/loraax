#include <fstream>
#include <vector>
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
// Default constructor
//
/******************************************************************************/
Airfoil::Airfoil ()
{
  _nb = 0;
  _n = 0;
  _xb.resize(0);
  _zb.resize(0);
  _x.resize(0);
  _z.resize(0);
  _s.resize(0);
  _xs.resize(0);
  _zs.resize(0);
  _y = 0.;
  _unit_transform = false;
  _cl = 0.;
  _cd = 0.;
  _cm = 0.;
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
  std::vector<std::string> xz;
  double x, z;
  int retval = 0;

  f.open(fname.c_str());
  if (not f.is_open()) { return 1; }

  // Start reading coordinates or skip if header line

  _nb = 0;
  _xb.resize(0);
  _zb.resize(0);
  std::getline(f, line); 
  xz = split_string(bracket_name(line));
  if ( (xz.size() == 2) && (string2double(xz[0], x) == 0) &&
       (string2double(xz[1], z) == 0) )
  {
    _nb += 1;
    _xb.push_back(x);
    _zb.push_back(z);
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
      _xb.push_back(x);
      _zb.push_back(z);
    } 
  }
    
  f.close();
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
  int n;
  unsigned int i;

  if (des.size() != 4)
    return 1;

  naca_4_digit(des.c_str(), &npointside, x, z, &n); 
  _nb = n;
  _xb.resize(_nb);
  _zb.resize(_nb);
  for ( i = 0; i < _nb; i++ )
  {
    _xb[i] = x[i];
    _zb[i] = z[i];
  }

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
  int n, stat;
  unsigned int i;

  if (des.size() != 5)
    return 1;

  naca_5_digit(des.c_str(), &npointside, x, z, &n, &stat); 
  if (stat != 0)
    return 2;

  _nb = n;
  _xb.resize(_nb);
  _zb.resize(_nb);
  for ( i = 0; i < _nb; i++ )
  {
    _xb[i] = x[i];
    _zb[i] = z[i];
  }

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
  if (x.size() != z.size())
    return 1;

  _nb = x.size();
  _xb = x;
  _zb = z;

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
  unsigned int i;
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
  _xb.resize(_nb);
  _zb.resize(_nb);
  foil1.smoothedCoordinates(x1, z1);
  foil2.smoothedCoordinates(x2, z2);
  for ( i = 0; i < _nb; i++ )
  {
    _xb[i] = (1.-interpfrac)*x1[i] + interpfrac*x2[i];
    _zb[i] = (1.-interpfrac)*z1[i] + interpfrac*z2[i];
  }

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
  unsigned int i;
  double dx1, dx2, dz1, dz2;
  std::vector<double> xrev(_nb), zrev(_nb);

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
    dx1 = _xb[i-1] - _xb[i];
    dx2 = _xb[i+1] - _xb[i];
    dz1 = _zb[i-1] - _zb[i];
    dz2 = _zb[i+1] - _zb[i];
    cross += (dz2*dx1 - dx2*dz1);
  }
  if (cross < 0.) { return; }

  // Reverse ordering for clockwise airfoil

  for ( i = 0; i < _nb; i++ )
  {
    xrev[i] = _xb[_nb-i-1];
    zrev[i] = _zb[_nb-i-1];
  }
  _xb = xrev;
  _zb = zrev;
}

/******************************************************************************/
//
// Fits a spline to buffer coordinates
//
/******************************************************************************/
int Airfoil::splineFit ()
{
  double x[_nb], z[_nb], s[_nb], xs[_nb], zs[_nb];
  unsigned int i;
  int n;

  if (_nb == 0)
  {
#ifdef DEBUG
    print_warning("Airfoil::splineFit", "Airfoil not yet loaded.");
#endif
    return 1;
  }

  n = _nb;
  _s.resize(_nb);
  _xs.resize(_nb);
  _zs.resize(_nb);
  for ( i = 0; i < _nb; i++ )
  {
    x[i] = _xb[i];
    z[i] = _zb[i];
  }
  xfoil_spline_coordinates(x, z, &n, s, xs, zs);

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
  double x[_nb], z[_nb], s[_nb], xs[_nb], zs[_nb];
  unsigned int i;
  int n;

  if (! splined())
  {
#ifdef DEBUG
    print_warning("Airfoil::splineInterp", "Must call splineFit first.");
#endif
    return 1;
  }

  n = _nb;
  for ( i = 0; i < _nb; i++ )
  {
    x[i] = _xb[i];
    z[i] = _zb[i];
    s[i] = _s[i];
    xs[i] = _xs[i];
    zs[i] = _zs[i];
  }
  xfoil_eval_spline(x, z, s, xs, zs, &n, &sc, &xc, &zc);

  return 0;
}

/******************************************************************************/
//
// Scales to unit chord and moves to origin
//
/******************************************************************************/
int Airfoil::unitTransform ()
{
  double xle, zle, xte, chord;
  double x[_nb], z[_nb], s[_nb], xs[_nb], zs[_nb];
  unsigned int i;
  int n;

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

  n = _nb;
  for ( i = 0; i < _nb; i++ )
  {
    x[i] = _xb[i];
    z[i] = _zb[i];
    s[i] = _s[i];
    xs[i] = _xs[i];
    zs[i] = _zs[i];
  }
  xfoil_lefind(x, z, s, xs, zs, &n, &_sle, &xle, &zle);

  xte = -1E+06;
  for ( i = 0; i < _nb; i++ )
  {
    if (_xb[i] > xte) { xte = _xb[i]; }
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
    _xb[i] -= xle;
    _zb[i] -= zle;
    _xb[i] /= chord;
    _zb[i] /= chord;
  }

  // Refit spline to transformed coordinates

  splineFit();
  for ( i = 0; i < _nb; i++ )
  {
    x[i] = _xb[i];
    z[i] = _zb[i];
    s[i] = _s[i];
    xs[i] = _xs[i];
    zs[i] = _zs[i];
  }
  xfoil_lefind(x, z, s, xs, zs, &n, &_sle, &xle, &zle);

  _unit_transform = true;

  return 0;
}

/******************************************************************************/
//
// Generates smoothed coordinates. Returns 1 and does nothing if buffer airfoil
// has not been loaded. Returns 2 and does nothing if unitTransform has not been
// called.
//
/******************************************************************************/
int Airfoil::smoothPaneling ( const xfoil_geom_options_type & geom_opts )
{
  unsigned int i;
  int npointin, npointout;
  double xin[_nb], zin[_nb];
  double xout[geom_opts.npan], zout[geom_opts.npan];

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

  for ( i = 0; i < _nb; i++ )
  {
    xin[i] = _xb[i]; 
    zin[i] = _zb[i]; 
  }

  npointin = _nb;
  npointout = geom_opts.npan;
  smooth_paneling(xin, zin, &npointin, &npointout, &geom_opts, xout, zout);

  _n = npointout;
  _x.resize(_n);
  _z.resize(_n);
  for ( i = 0; i < _n; i++ )
  {
    _x[i] = xout[i];
    _z[i] = zout[i];
  }

  return 0;
}

/******************************************************************************/
//
// Get airfoil data
//
/******************************************************************************/
unsigned int Airfoil::nBuffer () const { return _nb; }
unsigned int Airfoil::nSmoothed () const { return _n; }
void Airfoil::bufferCoordinates ( std::vector<double> & xb,
                                  std::vector<double> & zb ) const
{
#ifdef DEBUG
  if (_nb == 0)
  {
    print_warning("Airfoil::bufferCoordinates",
                  "Airfoil not yet loaded.");
  }
#endif

  xb = _xb;
  zb = _zb;
}

void Airfoil::smoothedCoordinates ( std::vector<double> & x,
                                    std::vector<double> & z ) const
{
#ifdef DEBUG
  if (_n == 0)
  {
    print_warning("Airfoil::smoothedCoordinates",
                  "Smoothed coordinates not yet available.");
  }
#endif

  x = _x;
  z = _z;
}
