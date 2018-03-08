#include <fstream>
#include <vector>
#include "util.h"
#include "airfoil.h"

/******************************************************************************/
//
// Airfoil class. Defines wing slice where BL equations are solved.
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
  _chord = 0.;
  _xb.resize(0);
  _zb.resize(0);
  _x.resize(0);
  _z.resize(0);
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
    if ( (xz.size() != 2) || (string2double(xz[0], x) == 0) ||
         (string2double(xz[1], z) == 0) )	// Read error
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

  return retval;
}

/******************************************************************************/
//
// Puts buffer coordinates in counterclockwise order. Does nothing if already
// in cc order.
//
/******************************************************************************/
void Airfoil::ccOrderCoordinates ()
{
  double cross;
  unsigned int i;
  double dx1, dx2, dz1, dz2;
  std::vector<double> xrev(_nb), zrev(_nb);

#ifndef NDEBUG
  if (_nb == 0)
  {
    conditional_stop(1, "Airfoil::ccOrderCoordinates",
                     "Airfoil not yet loaded.");
  }
#endif

  // Cross product e2 x e1 integral will be in -y direction for cc ordering

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
    std::cout << _nb-i-1 << std::endl;
    xrev[i] = _xb[_nb-i-1];
    zrev[i] = _zb[_nb-i-1];
  }
  _xb = xrev;
  _zb = zrev;
}
