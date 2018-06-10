#define _USE_MATH_DEFINES

#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "algorithms.h"
#include "util.h"
#include "settings.h"
#include "transformations.h"
#include "geometry.h"
#include "sectional_object.h"
#include "airfoil.h"
#include "section.h"
#include "vertex.h"
#include "quadpanel.h"
#include "tripanel.h"
#include "wake.h"
#include "wake_strip.h"
#include "wing.h"

// Data for optimizing spanwise spacing

std::vector<double> nom_stations, new_stations, s_wing;
std::vector<int> fixed_stations;

/******************************************************************************/
//
// Wing class. Contains sections, airfoils, panels, wake, etc.
// 
/******************************************************************************/

/******************************************************************************/
//
// Evaluates a possible combination of fixed stations and returns the sum of
// squared distance moved from the nominal spacing distribution
//
/******************************************************************************/
double evaluate_combination ( std::vector<int> & comb, unsigned int nfixed )
{
  double dist2;
  unsigned int i;

#ifdef DEBUG
  if (comb.size() != fixed_stations.size())
    conditional_stop(1, "evaluate_combination", "Vector size mismatch.");
#endif

  dist2 = 0;
  for ( i = 0; i < nfixed; i++ )
  {
    dist2 += std::pow(s_wing[i+1]-nom_stations[comb[i]], 2.);
  }

  return dist2;
}

/******************************************************************************/
//
// Evaluates all possible combinations of sections that could be fixed to user
// specified locations and saves the best one.
//
// Warning: this function may blow up the stack depending on the number of
// user fixed sections and spanwise stations, if not compiled with -O2 or -Os.
// If compiled with those optimizations, g++ should make this not really
// recursive "under the hood." This can be checked by compiling with the -S
// flag, opening with a text editor, and observing that this function does not
// call itself in the optimized code.
//
/******************************************************************************/
void optimize_fixed_stations ( std::vector<int> & comb, unsigned int nfixed,
                               unsigned int nsta, std::vector<int> & bestcomb,
                               double & bestval, unsigned int & evals )
{
  int i;
  double val;

  // Create initial combination

  if (comb.size() == 0)
  {
    for ( i = 0; i < int(nfixed); i++ )
    {
      comb.push_back(i+1);
    }
    comb[nfixed-1] = nfixed-1;
    bestcomb = comb;
    bestval = 1.E+06;
    evals = 0;
  }

  // Scan backwards through the list of possible combinations

  for ( i = nfixed-1; i >= 0; i-- )
  {
    if (comb[i] < int(nsta-2-(nfixed-(i+1))))
    {
      comb[i] += 1;
      if (i == int(nfixed-1))
      {
        val = evaluate_combination(comb, nfixed);
        evals += 1;
        if (val < bestval)
        {
          bestcomb = comb;
          bestval = val;
        }
      }
      optimize_fixed_stations(comb, nfixed, nsta, bestcomb, bestval, evals);
      // Return here lets g++ perform tail-recursion optimization
      return;
    }
    else
    {
      if (i > 0)
        comb[i] = comb[i-1] + 1;
    }
  }
} 

/******************************************************************************/
//
// Objective function to optimize spacing after setting fixed stations
//
/******************************************************************************/
double smooth_nonfixed_stations ( const std::vector<double> & ds_sta )
{
  unsigned int nsecs, nsta, i, j;
  double space, nom_space, objval;
  
  nsta = nom_stations.size();
  nsecs = s_wing.size();

  // Construct stations from deltas and fixed stations

  new_stations[0] = 0.;
  new_stations[nsta-1] = s_wing[nsecs-1];
  j = 0;
  for ( i = 1; i < nsta-1; i++ )
  {
    new_stations[i] = nom_stations[i] + ds_sta[i-1-j];
    if (j < fixed_stations.size())
    {
      if (int(i) == fixed_stations[j])
      {
        new_stations[i] = s_wing[j+1];
        j += 1;
      }
    }
  }

  // Penalize differences in spacing compared to nominal

  objval = 0.;
  for ( i = 1; i < nsta; i++ )
  {
    nom_space = nom_stations[i] - nom_stations[i-1];
    space = new_stations[i] - new_stations[i-1];
    if (space < 0.)
      objval += 5.*std::pow((space-nom_space)/nom_space, 2.);
    else
      objval += std::pow((space-nom_space)/nom_space, 2.);
  }

  return objval;
}

/******************************************************************************/
//
// Optimize spanwise spacing to put sections where the user asked while
// remaining as close as possible to the nominal distribution
//
/******************************************************************************/
std::vector<double> Wing::adjustSpacing (
                                const std::vector<double> & nom_stations ) const
{
  unsigned int i, j, nsta, nsecs;
  int nadjust, nfixed;
  std::vector<double> ds_sta0, ds_sta_opt;
  std::vector<int> comb;
  double span, fmin;
  unsigned int nsteps, fevals;
  conjgrad_options_type conjopt;
  double (*objfunc)(const std::vector<double> & x);

  // Initialize data

  nsecs = s_wing.size();
  span = s_wing[nsecs-1];
  nsta = nom_stations.size();

  // Optimize placement of fixed stations

  nfixed = nsecs-2;
  optimize_fixed_stations(comb, nfixed, nsta, fixed_stations, fmin, fevals);

  // Optimize placement of non-fixed stations

  conjopt.tol = 1.E-09;
  conjopt.h = 1e-12;
  conjopt.dxmax = span/double(nsta-1);
  conjopt.maxit = 2000;
#ifdef DEBUG
  conjopt.display_progress = true;
#else
  conjopt.display_progress = false;
#endif

  nadjust = nsta-2 - (nsecs-2);
  if (nadjust < 1)
    conditional_stop(1, "Wing::adjustSpacing", "Not enough spanwise points.");
  ds_sta0.resize(nadjust); 
  ds_sta_opt.resize(nadjust);
  new_stations.resize(nsta);
  objfunc = &smooth_nonfixed_stations;
  fevals = 0;
  for ( i = 0; int(i) < nadjust; i++ )
  {
    ds_sta0[i] = 0.;
  }
  conjgrad_search(ds_sta_opt, fmin, nsteps, fevals, objfunc, ds_sta0,
                  conjopt);

  // Populate new_stations vector

  new_stations[0] = nom_stations[0];
  new_stations[nsta-1] = nom_stations[nsta-1];
  j = 0;
  for ( i = 1; i < nsta-1; i++ )
  {
    new_stations[i] = nom_stations[i] + ds_sta_opt[i-1-j];
    if (j < fixed_stations.size())
    {
      if (int(i) == fixed_stations[j])
      {
        new_stations[i] = s_wing[j+1];
        j += 1;
      }
    }
  }

  return new_stations;
}

/******************************************************************************/
//
// Computes wing area and MAC from section data and writes to stdout. Also
// writes Re based on MAC.
//
/******************************************************************************/
void Wing::computeAreaMAC ( const std::vector<Section> &
                            sorted_user_sections ) const
{
  unsigned int i, nsecs;
  double span, area, cenx, dceny, ceny, cenz, cbar, cbarx;
  double xT1, xL1, y1, z1, xT2, xL2, y2, z2;
  double mT, mL, bT, bL, mp, mm, bp, bm;
  double ar, tr;

  nsecs = sorted_user_sections.size();
  span = 2.*sorted_user_sections[nsecs-1].y();
  tr = sorted_user_sections[nsecs-1].chord() / sorted_user_sections[0].chord();
  area = 0.;
  cenx = 0.;
  ceny = 0.;
  cenz = 0.;
  cbar = 0.;
  for ( i = 0; i < nsecs-1; i++ )
  {
    xL1 = sorted_user_sections[i].xle();
    xT1 = sorted_user_sections[i].xle() + sorted_user_sections[i].chord();
    y1 = sorted_user_sections[i].y();
    z1 = sorted_user_sections[i].zle();
    xL2 = sorted_user_sections[i+1].xle();
    xT2 = sorted_user_sections[i+1].xle() + sorted_user_sections[i+1].chord();
    y2 = sorted_user_sections[i+1].y();
    z2 = sorted_user_sections[i+1].zle();

    if (y2 - y1 < 1.E-12)
      continue;

    // Slopes and intercepts of LE and TE edges

    mT = (xT2 - xT1) / (y2 - y1);
    mL = (xL2 - xL1) / (y2 - y1);
    bT = xT1 - mT*y1;
    bL = xL1 - mL*y1;

    // Added and subtracted terms

    mp = mT + mL;
    mm = mT - mL;
    bp = bT + bL;
    bm = bT - bL;

    // Area, centroid, and MAC increments. Note that the increments to cenx,
    // ceny, and cbar here are weighted by area increment.

    area += 0.5*(mT - mL)*(std::pow(y2,2.) - std::pow(y1,2.))
         +  (bT - bL)*(y2 - y1);  
    cenx += 0.5 * (
            1./3.*mp*mm*(std::pow(y2,3.) - std::pow(y1,3.))
         +  0.5*(mp*bm + mm*bp)*(std::pow(y2,2.) - std::pow(y1,2.))
         +  bp*bm*(y2 - y1) );
    dceny = 1./3.*mm*(std::pow(y2,3.) - std::pow(y1,3.))
          +  0.5*bm*(std::pow(y2,2.) - std::pow(y1,2.));
    ceny += dceny;
    cenz += (dceny - y1) / (y2 - y1) * (z2 - z1) + z1;
    cbar += 1./3.*std::pow(mm,2.)*(std::pow(y2,3.) - std::pow(y1,3.))
         +  mm*bm*(std::pow(y2,2.) - std::pow(y1,2.))
         +  std::pow(bm,2.)*(y2 - y1);
  }
  if (area > 0.)
  {
    cenx /= area;
    ceny /= area;
    cenz /= area;
    cbar /= area;
    cbarx = cenx - 0.25*cbar;
    ar = std::pow(0.5*span,2.) / (0.5*area);
  }
  area *= 2.;	// Mirror image

  std::cout.setf(std::ios_base::scientific);
  std::cout << "Geometric information for wing " << _name << ":" << std::endl;
  std::cout << "  Span: "
            << std::setprecision(5) << span << std::endl;
  std::cout << "  Planform area: "
            << std::setprecision(5) << area << std::endl;
  std::cout << "  Mean aerodynamic chord: "
            << std::setprecision(5) << cbar << std::endl;
  std::cout << "  Halfspan MAC location: "
            << std::setprecision(5) << cbarx << ",  " << ceny << ", " 
            << cenz << std::endl;
  std::cout << "  Aspect ratio: "
            << std::setprecision(5) << ar << std::endl;
  std::cout << "  Taper ratio: "
            << std::setprecision(5) << tr << std::endl;
  std::cout << "  Re based on MAC: "
            << std::setprecision(5) << rhoinf * uinf * cbar / muinf
            << std::endl;
}

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
Wing::Wing ()
{
  _name = "";
  _nchord = 0;
  _nspan = 0;
  _ntipcap = 3;
  _lesprat = 1.;
  _tesprat = 1.;
  _rootsprat = 1.;
  _tipsprat = 1.;
  _sections.resize(0);
  _foils.resize(0);
  _verts.resize(0);
  _tipverts.resize(0);
  _quads.resize(0);
  _tris.resize(0);
  _panels.resize(0);
  _wakestrips.resize(0);
}

/******************************************************************************/
//
// Set/get name
//
/******************************************************************************/
void Wing::setName ( const std::string & name ) { _name = name; }
const std::string & Wing::name () const { return _name; }

/******************************************************************************/
//
// Set discretization and spacing info
//
/******************************************************************************/
void Wing::setDiscretization ( unsigned int nchord, unsigned int nspan,
                               const double & lesprat, const double & tesprat,
                               const double & rootsprat,
                               const double & tipsprat )
{
  _nchord = nchord;
  _nspan = nspan;
  _lesprat = lesprat;
  _tesprat = tesprat;
  _rootsprat = rootsprat;
  _tipsprat = tipsprat;
}

/******************************************************************************/
//
// Sets user airfoils from a list and sorts them in ascending y order.
//
/******************************************************************************/
void Wing::setAirfoils ( std::vector<Airfoil> & foils )
{
  SectionalObject *sorted_foils[foils.size()];
  unsigned int i, j, nfoils;

  nfoils = foils.size();
  if (nfoils == 1)
  {
    _foils = foils;
    return;
  }

  // Sort airfoils

  for ( i = 0; i < nfoils; i++ )
  {
    sorted_foils[i] = &foils[i];
  }
  sort_sections(sorted_foils, nfoils);
  for ( i = 0; i < nfoils; i++ )
  {
    for ( j = 0; j < nfoils; j++ )
    {
      if (foils[j].y() == sorted_foils[i]->y())
      {
        _foils.push_back(foils[j]);
        continue;
      }
    } 
  }

#ifdef DEBUG
  if (foils.size() < nfoils)
    conditional_stop(1, "Wing::setAirfoils", "Error sorting airfoils.");
#endif
}

/******************************************************************************/
//
// Set sections based on user section inputs and spacing. User sections define
// the planform shape, and the _sections variable defines the final
// discretization.
//
/******************************************************************************/
int Wing::setupSections ( std::vector<Section> & user_sections )
{
  SectionalObject *sorted_sections[user_sections.size()];
  std::vector<Section> sorted_user_sections;
  unsigned int i, j, nsecs, nfoils;
  Eigen::Vector2d secvec;
  double a4, a5, rootsp, tipsp, space, tipy, deltay, deltas, sfrac;
  double xle, zle, y, chord, twist, roll;
  std::vector<double> final_stations, foil_positions;

  // Sort sections

  nsecs = user_sections.size();
  for ( i = 0; i < nsecs; i++ )
  {
    sorted_sections[i] = &user_sections[i];
  }
  sort_sections(sorted_sections, nsecs);
  for ( i = 0; i < nsecs; i++ )
  {
    for ( j = 0; j < nsecs; j++ )
    {
      if (user_sections[j].y() == sorted_sections[i]->y())
      {
        sorted_user_sections.push_back(user_sections[j]);
        continue;
      }
    }
  }

#ifdef DEBUG
  if (sorted_user_sections.size() < nsecs)
  {
    conditional_stop(1, "Wing::setupSections", "Error sorting sections.");
    return 1;
  }
#endif

  s_wing.resize(nsecs);
  s_wing[0] = 0.;
  for ( i = 1; i < nsecs; i++ )
  {
    secvec[0] = sorted_user_sections[i].y() - sorted_user_sections[i-1].y();
    secvec[1] = sorted_user_sections[i].zle() - sorted_user_sections[i-1].zle();
    s_wing[i] = s_wing[i-1] + secvec.norm();
    if (s_wing[i] == s_wing[i-1])
    {
      conditional_stop(1, "Wing::setupSections",
                       "Wing has duplicate Y sections.");
      return 1;
    }
    sorted_user_sections[i].setRoll(atan(secvec[1]/secvec[0])*180./M_PI);
  }

  // Compute and print wing geometric info

  computeAreaMAC(sorted_user_sections);

  // Compute nominal discretized section locations (stations) from spacing

  rootsp = _rootsprat * s_wing[nsecs-1] / double(_nspan-1);
  tipsp = _tipsprat * s_wing[nsecs-1] / double(_nspan-1);
  opt_tanh_spacing(_nspan, s_wing[nsecs-1], rootsp, tipsp, a4, a5);
  nom_stations.resize(_nspan);
  nom_stations[0] = 0.; 
  for ( i = 1; i < _nspan; i++ )
  {
    space = tanh_spacing(i-1, a4, a5, _nspan, s_wing[nsecs-1], rootsp, tipsp);
    nom_stations[i] = nom_stations[i-1] + space;
  }

  // Optimize stations to be as close as possible to nominal but lining up with
  // user sections

  if (nsecs > 2)
    final_stations = adjustSpacing(nom_stations);
  else
    final_stations = nom_stations;

  // Add airfoils as needed to ensure the entire span is covered

  nfoils = _foils.size();
  tipy = sorted_user_sections[nsecs-1].y();
  if (nfoils < 1)
  {
    conditional_stop(1, "Wing::setupSections",
                     "At least one airfoil is required.");
    return 1;
  }
  else if (nfoils == 1)
  {
    // If a single airfoil is provided, modify _foils vector so it has
    // identical airfoils at the root and tip

    _foils.push_back(_foils[0]);
    _foils[0].setY(0.);
    _foils[1].setY(tipy);
    nfoils += 1;
  } 
  else
  {
    // Add foils at the root and tip if none were specified there

    if (_foils[0].y() > 0)
    {
      _foils.insert(_foils.begin(), _foils[0]);
      nfoils += 1;
      _foils[0].setY(0.);
    }
    if (_foils[nfoils-1].y() < tipy)
    {
      _foils.push_back(_foils[nfoils-1]);
      nfoils += 1;
      _foils[nfoils-1].setY(tipy);
    }
    else if (_foils[nfoils-1].y() > tipy)
    {
      conditional_stop(1, "Wing::setupSections",
                       "Airfoil is placed beyond the wingtip.");
      return 1;
    }
  } 

  // Determine s_wing spanwise positions for airfoils

  foil_positions.resize(nfoils);
  for ( i = 0; i < nfoils; i++ )
  {
    for ( j = 0; j < nsecs-1; j++ )
    {
      if ( (_foils[i].y() >= sorted_user_sections[j].y()) &&
           (_foils[i].y() <= sorted_user_sections[j+1].y()) )
      {
        deltay = sorted_user_sections[j+1].y() - sorted_user_sections[j].y(); 
        sfrac = (_foils[i].y() - sorted_user_sections[j].y()) / deltay;
        deltas = s_wing[j+1] - s_wing[j];
        foil_positions[i] = s_wing[j] + sfrac*deltas;
        break;
      }
    }
  }

  // Create sections

  _sections.resize(_nspan);
  for ( i = 0; i < _nspan; i++ )
  {
    // Set airfoil coordinates

    for ( j = 0; j < nfoils-1; j++ )
    {
      if (std::abs(final_stations[i]-foil_positions[j]) < 1e-12)
      {
        _sections[i].airfoil().interpCoordinates(_foils[j], _foils[j], 0.);
        break;
      }
      else if ( (final_stations[i] > foil_positions[j]) &&
                (final_stations[i] < foil_positions[j+1]) )
      {
        sfrac = (final_stations[i] - foil_positions[j]) / 
                (foil_positions[j+1] - foil_positions[j]);
        _sections[i].airfoil().interpCoordinates(_foils[j], _foils[j+1], sfrac);  
        break;
      }
      else if ( (j == nfoils-2) &&
                (std::abs(final_stations[i]-foil_positions[j+1]) < 1e-12) )
      {
        _sections[i].airfoil().interpCoordinates(_foils[j+1], _foils[j+1], 0.);
        break;
      }
      else if (j == nfoils-2)
      {
        conditional_stop(1, "Wing::setupSections",
                         "Could not find interpolant airfoils.");
        return 1;
      }
    }
    _sections[i].airfoil().ccwOrderCoordinates();
    _sections[i].airfoil().splineFit();
    _sections[i].airfoil().unitTransform();
    _sections[i].airfoil().setXfoilOptions(xfoil_run_opts, xfoil_geom_opts);
    _sections[i].airfoil().smoothPaneling();

    // Set section position, orientation, and scale

    for ( j = 0; j < nsecs-1; j++ )
    {
      if (std::abs(final_stations[i]-s_wing[j]) < 1e-12)
      {
        xle = sorted_user_sections[j].xle();
        zle = sorted_user_sections[j].zle();
        y = sorted_user_sections[j].y();
        chord = sorted_user_sections[j].chord();
        twist = sorted_user_sections[j].twist();
        roll = sorted_user_sections[j].roll();
        break;
      }
      else if ( (final_stations[i] > s_wing[j]) &&
                (final_stations[i] < s_wing[j+1]) )
      {
        sfrac = (final_stations[i] - s_wing[j]) /
                (s_wing[j+1] - s_wing[j]);
          xle = (1.-sfrac)*sorted_user_sections[j].xle() +
                    sfrac*sorted_user_sections[j+1].xle();
            y = (1.-sfrac)*sorted_user_sections[j].y() +
                    sfrac*sorted_user_sections[j+1].y();
          zle = (1.-sfrac)*sorted_user_sections[j].zle() +
                    sfrac*sorted_user_sections[j+1].zle();
        chord = (1.-sfrac)*sorted_user_sections[j].chord() +
                    sfrac*sorted_user_sections[j+1].chord();
        twist = (1.-sfrac)*sorted_user_sections[j].twist() +
                    sfrac*sorted_user_sections[j+1].twist();
         roll = sorted_user_sections[j+1].roll();
        break;
      }
      else if ( (j == nsecs-2) &&
                (std::abs(final_stations[i]-s_wing[j+1]) < 1e-12) )
      {
        xle = sorted_user_sections[j+1].xle();
        zle = sorted_user_sections[j+1].zle();
        y = sorted_user_sections[j+1].y();
        chord = sorted_user_sections[j+1].chord();
        twist = sorted_user_sections[j+1].twist();
        roll = sorted_user_sections[j+1].roll();
        break;
      }
      else if (j == nsecs-2)
      {
        conditional_stop(1, "Wing::setupSections",
                         "Could not find interpolant sections.");
        return 1;
      }
    }
    _sections[i].setGeometry(xle, y, zle, chord, twist, roll);

    // Set vertices from spacing distribution

    _sections[i].setVertices(_nchord, _lesprat, _tesprat);
    
    // Set Reynolds number and Mach number if viscous
    
    if (viscous) 
    {
      _sections[i].computeReynoldsNumber(rhoinf, uinf, muinf);
      if (compressible)
        _sections[i].setMachNumber(minf);
      else
        _sections[i].setMachNumber(0.0);
    }
  }

  return 0; 
}

/******************************************************************************/
//
// Creates panels and surface vertex pointers. Increments next_global_vertidx
// and next_global_elemidx
//
/******************************************************************************/
void Wing::createPanels ( int & next_global_vertidx, int & next_global_elemidx )
{
  unsigned int i, j, vcounter, qcounter, tcounter, ntri, nquad, right;
  double phin, phi, r;
  Eigen::Matrix3d trans, T1;
  Eigen::Vector3d cen, r0, rb, ri, point, norm, tang, tangb;
  Eigen::Vector3d tanl, tanr, tanf, tanb, tan;
  double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;

  // Set vertex pointers on top and bottom surfaces

  _verts.resize(_nspan*(2*_nchord-1) + (_ntipcap-2)*(_nchord-2));
  vcounter = 0;
  for ( i = 0; i < _nspan; i++ )
  {
    for ( j = 0; j < 2*_nchord-1; j++ )
    {
      _sections[i].vert(j).setIdx(next_global_vertidx);
      _verts[vcounter] = &_sections[i].vert(j);
      vcounter += 1;
      next_global_vertidx += 1;
    }
  }

  // Determine number of tris and quads

  ntri = 2*(_ntipcap-1);
  nquad = (_nspan-1)*(2*_nchord-2) + (_ntipcap-1)*(_nchord-3);
 
  // Create quad panels on top/bottom surfaces

  _quads.resize(nquad);
  qcounter = 0;
  _panels.resize(_nspan-1 + (_ntipcap-1)/2);
  for ( i = 0; i < _nspan-1; i++ )
  {
    _panels[i].resize(2*_nchord-2);
    for ( j = 0; j < 2*_nchord-2; j++ )
    {
      _quads[qcounter].setIdx(next_global_elemidx);
      _quads[qcounter].addVertex(&_sections[i].vert(j));
      _quads[qcounter].addVertex(&_sections[i+1].vert(j));
      _quads[qcounter].addVertex(&_sections[i+1].vert(j+1));
      _quads[qcounter].addVertex(&_sections[i].vert(j+1));

	  // Compute surface tangent vector at centroid

	  tanl << _sections[i].vert(j+1).x() - _sections[i].vert(j).x(),
	          _sections[i].vert(j+1).y() - _sections[i].vert(j).y(),
	          _sections[i].vert(j+1).z() - _sections[i].vert(j).z();
	  tanr << _sections[i+1].vert(j+1).x() - _sections[i+1].vert(j).x(),
	          _sections[i+1].vert(j+1).y() - _sections[i+1].vert(j).y(),
	          _sections[i+1].vert(j+1).z() - _sections[i+1].vert(j).z();
	  tan = 0.5*(tanl + tanr);
	  tan /= tan.norm();
	  _quads[qcounter].setTangent(tan);

      _panels[i][j] = &_quads[qcounter];
      qcounter += 1;
      next_global_elemidx += 1;
    }
  }

  // Create tip cap vertices. Currently, there's a lot of stuff in here that
  // does nothing useful, because I switched back to flat tips. The code is
  // kept here in case I want to go back to revolved tip caps again.

  _tipverts.resize(_ntipcap-2);
  for ( i = 1; i < _ntipcap-1; i++ )
  {
    _tipverts[i-1].resize(_nchord-2);
    for ( j = 1; j < _nchord-1; j++ )
    {
      // Arc angle

      phi = double(i)/double(_ntipcap-1)*180.;

      // Get normal vector from last spanwise section

      x1 = _sections[_nspan-2].vert(j).x();
      y1 = _sections[_nspan-2].vert(j).y();
      z1 = _sections[_nspan-2].vert(j).z();
      x2 = _sections[_nspan-1].vert(j).x();
      y2 = _sections[_nspan-1].vert(j).y();
      z2 = _sections[_nspan-1].vert(j).z();
      x3 = _sections[_nspan-1].vert(2*_nchord-2-j).x();
      y3 = _sections[_nspan-1].vert(2*_nchord-2-j).y();
      z3 = _sections[_nspan-1].vert(2*_nchord-2-j).z();
      x4 = _sections[_nspan-2].vert(2*_nchord-2-j).x();
      y4 = _sections[_nspan-2].vert(2*_nchord-2-j).y();
      z4 = _sections[_nspan-2].vert(2*_nchord-2-j).z();
      norm = quad_normal(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);

      // Transform from normal vector

      trans = transform_from_normal(norm(0), norm(1), norm(2));

      // Tangential vector is the projection of the tip section's normal vector
      // on this plane

      tang(0) = 0.;
      tang(1) = cos(_sections[_nspan-1].roll()*M_PI/180.);
      tang(2) = sin(_sections[_nspan-1].roll()*M_PI/180.);

      // Account for local swept dihedral angle

      tangb = trans*tang;
      phin = atan(tangb(0)/tangb(1))*180./M_PI; 
      T1 = euler_rotation(0., 0., -phin);

      // Center of revolution and radial vector in x-y plane

      cen(0) = 0.5*(_sections[_nspan-1].vert(j).x() + 
                    _sections[_nspan-1].vert(2*_nchord-2-j).x());
      cen(1) = 0.5*(_sections[_nspan-1].vert(j).y() + 
                    _sections[_nspan-1].vert(2*_nchord-2-j).y());
      cen(2) = 0.5*(_sections[_nspan-1].vert(j).z() + 
                    _sections[_nspan-1].vert(2*_nchord-2-j).z());
      r0(0) = _sections[_nspan-1].vert(j).x() - cen(0);
      r0(1) = _sections[_nspan-1].vert(j).y() - cen(1);
      r0(2) = _sections[_nspan-1].vert(j).z() - cen(2);
      /* Use this to enable rounded tip caps. Also would need to make _ntipcap
         an input again in that case.
      r = r0.norm();
      */
      r = 0.;

      // Radial vector in y-z plane

      rb << r*cos(phi*M_PI/180.), r*sin(phi*M_PI/180.), 0.;

      // Compute vertex location and create vertex

      ri = (T1*trans).transpose()*rb;
      point = cen + ri;
      _tipverts[i-1][j-1].setIdx(next_global_vertidx);
      _tipverts[i-1][j-1].setCoordinates(point(0), point(1), point(2));
      _verts[vcounter] = &_tipverts[i-1][j-1];
      vcounter += 1;
      next_global_vertidx += 1;
    }
  } 

  // Create tip panels wrapping around the tip from top TE to bottom TE.
  // First loop connects to existing vertices on last section. Note: don't
  // "connect" to vertices on top and bottom surfaces, because we don't want
  // the tip panels to be included in averaging to those vertices.

  _tris.resize(ntri);
  tcounter = 0;
  _panels[_nspan-1].resize(2*_nchord-2);
  for ( j = 0; j < 2*_nchord-2; j++ )
  {
    // Tri panels at TE and LE

    if (j == 0)
    {
      _tris[tcounter].setIdx(next_global_elemidx);
      _tris[tcounter].addVertex(&_sections[_nspan-1].vert(0),false);
      _tris[tcounter].addVertex(&_tipverts[0][0]);
      _tris[tcounter].addVertex(&_sections[_nspan-1].vert(1),false);
      _panels[_nspan-1][j] = &_tris[tcounter];
      tcounter += 1;
      next_global_elemidx += 1;
    }
    else if (j == _nchord-2)
    {
      _tris[tcounter].setIdx(next_global_elemidx);
      _tris[tcounter].addVertex(&_tipverts[0][_nchord-3]);
      _tris[tcounter].addVertex(&_sections[_nspan-1].vert(_nchord-1),false);
      _tris[tcounter].addVertex(&_sections[_nspan-1].vert(_nchord-2),false);
      _panels[_nspan-1][j] = &_tris[tcounter];
      tcounter += 1;
      next_global_elemidx += 1;
    }
    else if (j == _nchord-1)
    {
      _tris[tcounter].setIdx(next_global_elemidx);
      _tris[tcounter].addVertex(&_sections[_nspan-1].vert(_nchord-1),false);
      _tris[tcounter].addVertex(&_tipverts[_ntipcap-3][_nchord-3]);
      _tris[tcounter].addVertex(&_sections[_nspan-1].vert(_nchord),false);
      _panels[_nspan-1][j] = &_tris[tcounter];
      tcounter += 1;
      next_global_elemidx += 1;
    }
    else if (j == 2*_nchord-3)
    {
      _tris[tcounter].setIdx(next_global_elemidx);
      _tris[tcounter].addVertex(&_tipverts[_ntipcap-3][0]);
      _tris[tcounter].addVertex(&_sections[_nspan-1].vert(2*_nchord-2),false);
      _tris[tcounter].addVertex(&_sections[_nspan-1].vert(2*_nchord-3),false);
      _panels[_nspan-1][j] = &_tris[tcounter];
      tcounter += 1;
      next_global_elemidx += 1;
    }

    // Quad panels in between

    else if (j < _nchord-2)
    {
      _quads[qcounter].setIdx(next_global_elemidx);
      _quads[qcounter].addVertex(&_tipverts[0][j-1]);
      _quads[qcounter].addVertex(&_tipverts[0][j]);
      _quads[qcounter].addVertex(&_sections[_nspan-1].vert(j+1),false);
      _quads[qcounter].addVertex(&_sections[_nspan-1].vert(j),false);
      _panels[_nspan-1][j] = &_quads[qcounter];
      qcounter += 1;
      next_global_elemidx += 1;
    }
    else
    {
      _quads[qcounter].setIdx(next_global_elemidx);
      _quads[qcounter].addVertex(&_tipverts[_ntipcap-3][2*_nchord-3-j]);
      _quads[qcounter].addVertex(&_tipverts[_ntipcap-3][2*_nchord-3-j-1]);
      _quads[qcounter].addVertex(&_sections[_nspan-1].vert(j+1),false);
      _quads[qcounter].addVertex(&_sections[_nspan-1].vert(j),false);
      _panels[_nspan-1][j] = &_quads[qcounter];
      qcounter += 1;
      next_global_elemidx += 1;
    }
  }

  // Next layers involve only _tipverts except at TE and LE points.
  // As before, don't connect to top/bottom surfaces vertices.

  for ( i = 1; i < (_ntipcap-1)/2; i++ )
  {
    _panels[_nspan-1+i].resize(2*_nchord-2);
    for ( j = 0; j < 2*_nchord-2; j++ )
    {
      // Tri panels at TE and LE

      if (j == 0)
      {
        _tris[tcounter].setIdx(next_global_elemidx);
        _tris[tcounter].addVertex(&_sections[_nspan-1].vert(0),false);
        _tris[tcounter].addVertex(&_tipverts[i][0]);
        _tris[tcounter].addVertex(&_tipverts[i-1][0]);
        _panels[_nspan-1+i][j] = &_tris[tcounter];
        tcounter += 1;
        next_global_elemidx += 1;
      }
      else if (j == _nchord-2)
      {
        _tris[tcounter].setIdx(next_global_elemidx);
        _tris[tcounter].addVertex(&_tipverts[i][_nchord-3]);
        _tris[tcounter].addVertex(&_sections[_nspan-1].vert(_nchord-1),false);
        _tris[tcounter].addVertex(&_tipverts[i-1][_nchord-3]);
        _panels[_nspan-1+i][j] = &_tris[tcounter];
        tcounter += 1;
        next_global_elemidx += 1;
      }
      else if (j == _nchord-1)
      {
        _tris[tcounter].setIdx(next_global_elemidx);
        _tris[tcounter].addVertex(&_sections[_nspan-1].vert(_nchord-1),false);
        _tris[tcounter].addVertex(&_tipverts[_ntipcap-3-i][_nchord-3]);
        _tris[tcounter].addVertex(&_tipverts[_ntipcap-3-i+1][_nchord-3]);
        _panels[_nspan-1+i][j] = &_tris[tcounter];
        tcounter += 1;
        next_global_elemidx += 1;
      }
      else if (j == 2*_nchord-3)
      {
        _tris[tcounter].setIdx(next_global_elemidx);
        _tris[tcounter].addVertex(&_tipverts[_ntipcap-3-i][0]);
        _tris[tcounter].addVertex(&_sections[_nspan-1].vert(2*_nchord-2),false);
        _tris[tcounter].addVertex(&_tipverts[_ntipcap-3-i+1][0]);
        _panels[_nspan-1+i][j] = &_tris[tcounter];
        tcounter += 1;
        next_global_elemidx += 1;
      }

      // Quad panels in between

      else if (j < _nchord-2)
      {
        _quads[qcounter].setIdx(next_global_elemidx);
        _quads[qcounter].addVertex(&_tipverts[i][j-1]);
        _quads[qcounter].addVertex(&_tipverts[i][j]);
        _quads[qcounter].addVertex(&_tipverts[i-1][j]);
        _quads[qcounter].addVertex(&_tipverts[i-1][j-1]);
        _panels[_nspan-1+i][j] = &_quads[qcounter];
        qcounter += 1;
        next_global_elemidx += 1;
      }
      else
      {
        _quads[qcounter].setIdx(next_global_elemidx);
        _quads[qcounter].addVertex(&_tipverts[_ntipcap-3-i][2*_nchord-3-j]);
        _quads[qcounter].addVertex(&_tipverts[_ntipcap-3-i][2*_nchord-3-j-1]);
        _quads[qcounter].addVertex(&_tipverts[_ntipcap-3-i+1][2*_nchord-3-j-1]);
        _quads[qcounter].addVertex(&_tipverts[_ntipcap-3-i+1][2*_nchord-3-j]);
        _panels[_nspan-1+i][j] = &_quads[qcounter];
        qcounter += 1;
        next_global_elemidx += 1;
      }
    }
  }

  // Surface tangent vectors on tip panels

  for ( i = 0; i < (_ntipcap-1)/2; i++ )
  {
    for ( j = 0; j < 2*_nchord-2; j++ )
	{
      // Tri panels at TE and LE

      if ( (j == 0) || (j == _nchord-1) )
      {
		tanf = _panels[_nspan-1+i][j+1]->centroid()
			 - _panels[_nspan-1+i][j]->centroid();
		tan = tanf / tanf.norm();
      }
      else if ( (j == _nchord-2) || (j == 2*_nchord-3) )
      {
		tanb = _panels[_nspan-1+i][j]->centroid()
			 - _panels[_nspan-1+i][j-1]->centroid();
		tan = tanb / tanb.norm();
      }

      // Quad panels in between

      else if (j < _nchord-2)
      {
		tanf = _panels[_nspan-1+i][j+1]->centroid()
			 - _panels[_nspan-1+i][j]->centroid();
		tanb = _panels[_nspan-1+i][j]->centroid()
			 - _panels[_nspan-1+i][j-1]->centroid();
		tan = 0.5*(tanf + tanb);
		tan /= tan.norm();
      }

      _panels[_nspan-1+i][j]->setTangent(tan);
    }
  }
	  

  // Set panel neighbors (top and bottom surfaces only for now)
  // We don't add panel neighbors from top/bottom to tip caps, because there can
  // be very large changes in sizing across that boundary, which would produce
  // error in gradient calculations.

  for ( i = 0; i < _nspan-1; i++ )
  {
    for ( j = 0; j < 2*_nchord-2; j++ )
    {
      if (i > 0)
        _panels[i][j]->setLeftNeighbor(_panels[i-1][j]);
      if (i < _nspan-2)
        _panels[i][j]->setRightNeighbor(_panels[i+1][j]);
      if (j > 0)
        _panels[i][j]->setBackNeighbor(_panels[i][j-1]);
      if (j < 2*_nchord-3)
        _panels[i][j]->setFrontNeighbor(_panels[i][j+1]);
    }
  }

  // Neighbor panels on tip caps

  for ( i = 0; i < (_ntipcap-1)/2; i++ )
  {
    for ( j = 0; j < 2*_nchord-2; j++ )
    {
      if (i > 0)
        _panels[_nspan-1+i][j]->setLeftNeighbor(_panels[_nspan-1+i-1][j]);
      if (i < (_ntipcap-1)/2-1)
        _panels[_nspan-1+i][j]->setRightNeighbor(_panels[_nspan-1+i+1][j]);
      if ( (j > 0) && (j != _nchord-1) )
        _panels[_nspan-1+i][j]->setBackNeighbor(_panels[_nspan-1+i][j-1]);
      if ( (j < 2*_nchord-3) && (j != _nchord-2) )
        _panels[_nspan-1+i][j]->setFrontNeighbor(_panels[_nspan-1+i][j+1]);

      // Add right neighbor across the tip cut

      if (i == (_ntipcap-1)/2-1)
      {
        right = 2*_nchord-3-j;
        _panels[_nspan-1+i][j]->setRightNeighbor(_panels[_nspan-1+i][right]);
      }
    }
  }

  // Compute grid metrics

#pragma omp parallel for private(i,j)
  for ( i = 0; i < _nspan-1+(_ntipcap-1)/2; i++ )
  {
    for ( j = 0; j < 2*_nchord-2; j++ )
    {
      _panels[i][j]->computeGridTransformation();
    }
  }
}

/******************************************************************************/
//
// Sets up wake
//
/******************************************************************************/
void Wing::setupWake ( int & next_global_vertidx, int & next_global_elemidx )
{
  unsigned int i, j, nstream;
  std::vector<Vertex *> topteverts, botteverts;

  // Create vertices along trailing edge

  topteverts.resize(_nspan);
  botteverts.resize(_nspan);
  for ( i = 0; i < _nspan; i++ )
  {
    topteverts[i] = &_sections[i].vert(0); 
    botteverts[i] = &_sections[i].vert(2*_nchord-2); 
  }

  // Initialize wake

  _wake.initialize(topteverts, botteverts, next_global_vertidx,
                   next_global_elemidx);

  // Create wake strips

  nstream = wakesteps+1;
  _wakestrips.resize(_nspan-1);
  for ( i = 0; i < _nspan-1; i++ )
  {
    _wakestrips[i].setNPanels((nstream-1)*2+1);
    _wakestrips[i].setTEPanels(&_quads[i*(2*_nchord-2)],
                               &_quads[(i+1)*(2*_nchord-2)-1]);
    for ( j = 0; j < nstream-1; j++ )
    {
      _wakestrips[i].setPanelPointer(j*2, _wake.triPanel(i*(nstream-1)*2+j*2));
      _wakestrips[i].setPanelPointer(j*2+1,
                                         _wake.triPanel(i*(nstream-1)*2+j*2+1));
    }
    _wakestrips[i].setPanelPointer((nstream-1)*2, _wake.quadPanel(i));
  }
}

/******************************************************************************/
//
// Computes velocities and pressures on surface panels, and interpolates to
// vertices
//
/******************************************************************************/
void Wing::computeSurfaceQuantities ()
{
  unsigned int i, j, nverts;
  double s12, s1, s2;
  Eigen::Matrix3d A;
  Eigen::Vector3d x, b;
  Eigen::PartialPivLU<Eigen::Matrix3d> lu;
  Vertex * v0, * v1, * v2;

#pragma omp parallel for private(i,j)
  for ( i = 0; i < _nspan-1+(_ntipcap-1)/2; i++ )
  {
    for ( j = 0; j < 2*_nchord-2; j++ )
    {
      _panels[i][j]->computeVelocity(uinfvec);
      _panels[i][j]->computePressure(uinf, rhoinf, pinf, minf, compressible);
    }
  }

  // Interpolate to vertices

  nverts = _verts.size();
#pragma omp parallel for private(i)
  for ( i = 0; i < nverts; i++ )
  {
    _verts[i]->interpFromPanels();
  }

  /* Extrapolate to vertices at edges using quadratic fit. Note that the
     originally calculated values at edge vertices actually apply to the
     midpoint of the boundary between adjacent panels, since it is averaged only
     from these two panels. */

  // Top trailing edge
#pragma omp parallel for private(i,v0,v1,v2,s1,s12,s2,A,lu,j,b,x)
  for ( i = 0; i < _nspan-1; i++ )
  {
    v0 = &_sections[i].vert(0);
    v1 = &_sections[i].vert(1);
    v2 = &_sections[i].vert(2);
    s1 = v0->distance(*v1);
    s12 = 0.5*s1;
    s2 = s1 + v1->distance(*v2);
    A << std::pow(s12,2.), s12, 1.,
         std::pow(s1,2.),  s1,  1.,
         std::pow(s2,2.),  s2,  1.;
    lu.compute(A);

    for ( j = 0; j < Vertex::dataSize; j++ )
    { 
      b << v0->data(j), v1->data(j), v2->data(j);
      x = lu.solve(b);
      v0->setData(j, x(2));
    }
  }

  // Bottom trailing edge
#pragma omp parallel for private(i,v0,v1,v2,s1,s12,s2,A,lu,j,b,x)
  for ( i = 0; i < _nspan-1; i++ )
  {
    v0 = &_sections[i].vert(2*_nchord-2);
    v1 = &_sections[i].vert(2*_nchord-3);
    v2 = &_sections[i].vert(2*_nchord-4);
    s1 = v0->distance(*v1);
    s12 = 0.5*s1;
    s2 = s1 + v1->distance(*v2);
    A << std::pow(s12,2.), s12, 1.,
         std::pow(s1,2.),  s1,  1.,
         std::pow(s2,2.),  s2,  1.;
    lu.compute(A);

    for ( j = 0; j < Vertex::dataSize; j++ )
    { 
      b << v0->data(j), v1->data(j), v2->data(j);
      x = lu.solve(b);
      v0->setData(j, x(2));
    }
  }

  // Centerline
#pragma omp parallel for private(i,v0,v1,v2,s1,s12,s2,A,lu,j,b,x)
  for ( i = 1; i < 2*_nchord-2; i++ )
  {
    v0 = &_sections[0].vert(i);
    v1 = &_sections[1].vert(i);
    v2 = &_sections[2].vert(i);
    s1 = v0->distance(*v1);
    s12 = 0.5*s1;
    s2 = s1 + v1->distance(*v2);
    A << std::pow(s12,2.), s12, 1.,
         std::pow(s1,2.),  s1,  1.,
         std::pow(s2,2.),  s2,  1.;
    lu.compute(A);

    for ( j = 0; j < Vertex::dataSize; j++ )
    { 
      b << v0->data(j), v1->data(j), v2->data(j);
      x = lu.solve(b);
      v0->setData(j, x(2));
    }
  }

  // Tip
#pragma omp parallel for private(i,v0,v1,v2,s1,s12,s2,A,lu,j,b,x)
  for ( i = 1; i < 2*_nchord-2; i++ )
  {
    v0 = &_sections[_nspan-1].vert(i);
    v1 = &_sections[_nspan-2].vert(i);
    v2 = &_sections[_nspan-3].vert(i);
    s1 = v0->distance(*v1);
    s12 = 0.5*s1;
    s2 = s1 + v1->distance(*v2);
    A << std::pow(s12,2.), s12, 1.,
         std::pow(s1,2.),  s1,  1.,
         std::pow(s2,2.),  s2,  1.;
    lu.compute(A);

    for ( j = 0; j < Vertex::dataSize; j++ )
    { 
      b << v0->data(j), v1->data(j), v2->data(j);
      x = lu.solve(b);
      v0->setData(j, x(2));
    }
  }
}

/******************************************************************************/
//
// Access to verts and panels
//
/******************************************************************************/
unsigned int Wing::nVerts () const { return _verts.size(); }
unsigned int Wing::nQuads () const { return _quads.size(); }
unsigned int Wing::nTris () const { return _tris.size(); }
Vertex * Wing::vert ( unsigned int vidx )
{
#ifdef DEBUG
  if (vidx >= _verts.size())
    conditional_stop(1, "Wing::vert", "Index out of range.");
#endif

  return _verts[vidx];
}

QuadPanel * Wing::quadPanel ( unsigned int qidx )
{
#ifdef DEBUG
  if (qidx >= _quads.size())
    conditional_stop(1, "Wing::quadPanel", "Index out of range.");
#endif

  return &_quads[qidx];
}

TriPanel * Wing::triPanel ( unsigned int tidx )
{
#ifdef DEBUG
  if (tidx >= _tris.size())
    conditional_stop(1, "Wing::triPanel", "Index out of range.");
#endif

  return &_tris[tidx];
}

/******************************************************************************/
//
// Access to wake and wake strips
//
/******************************************************************************/
Wake & Wing::wake () { return _wake; }
unsigned int Wing::nWStrips () const { return _wakestrips.size(); }
WakeStrip * Wing::wStrip ( unsigned int wsidx )
{
#ifdef DEBUG
  if (wsidx >= _wakestrips.size())
    conditional_stop(1, "Wing::wStrip", "Index out of range.");
#endif

  return &_wakestrips[wsidx];
}

/******************************************************************************/
//
// Computes viscous forces (and skin friction, etc.) using Xfoil at sections
//
/******************************************************************************/
void Wing::computeBL ()
{
  unsigned int i, j, k;
  double weighttop, weightbot, var;

#pragma omp parallel for private(i)
  for ( i = 0; i < _nspan; i++ )
  {
    _sections[i].computeBL(uinfvec, rhoinf, alpha);
    if (not _sections[i].blConverged())
    {
      std::cout << "    Warning: Xfoil BL calculations did not converge "
                << "for section " << i+1 << "." << std::endl;
    }
  }

  // Inteprolate BL quantities to tip vertices
  
  for ( i = 1; i < _ntipcap-1; i++ )
  {
    weightbot = double(i) / double(_ntipcap-1);
    weighttop = 1. - weightbot;
    for ( j = 1; j < _nchord-1; j++ )
    {
      for ( k = Vertex::firstBLData; k < Vertex::dataSize; k++ )
      {
        var = weighttop*_sections[_nspan-1].vert(j).data(k)
            + weightbot*_sections[_nspan-1].vert(2*_nchord-2-j).data(k);
        _tipverts[i-1][j-1].setData(k, var);
      }
    }
  }
}

/******************************************************************************/
//
// Compute or access forces and moments. Compute part includes section forces
// and moments.
//
/******************************************************************************/
void Wing::computeForceMoment ( const double & sref, const double & lref,
                                const Eigen::Vector3d & momcen )
{
	unsigned int i, j;
	Eigen::Vector3d dfp, dfv, dmp, dmv, fp, fv, mp, mv;
	double qinf;
	
	fp << 0., 0., 0.;
	fv << 0., 0., 0.;
	mp << 0., 0., 0.;
	mv << 0., 0., 0.;
	
	// Compute on top and bottom surfaces
	
	for ( i = 0; i < _nspan-1+(_ntipcap-1)/2; i++ )
	{
		for ( j = 0; j < 2*_nchord-2; j++ )
		{ 
		_panels[i][j]->computeForceMoment(uinf, rhoinf, pinf, momcen, viscous,
		                                  dfp, dfv, dmp, dmv);
		fp += dfp;
		fv += dfv;
		mp += dmp;
		mv += dmv;
		}
	}
	
	// Account for mirror image
	
	fp(0)  *= 2.;
	fp(1)  =  0.;
	fp(2)  *= 2.;
	fv(0)  *= 2.;
	fv(1)  =  0.;
	fv(2)  *= 2.;
	mp(0) =  0.;
	mp(1) *= 2.;
	mp(2) =  0.;
	mv(0) =  0.;
	mv(1) *= 2.;
	mv(2) =  0.;
	
	// Convert to wind frame
	
	qinf = 0.5*rhoinf*std::pow(uinf, 2.);
	_liftp = -fp(0)*sin(alpha*M_PI/180.) + fp(2)*cos(alpha*M_PI/180.);
	_liftv = -fv(0)*sin(alpha*M_PI/180.) + fv(2)*cos(alpha*M_PI/180.);
	_dragp =  fp(0)*cos(alpha*M_PI/180.) + fp(2)*sin(alpha*M_PI/180.);
	_dragv =  fv(0)*cos(alpha*M_PI/180.) + fv(2)*sin(alpha*M_PI/180.);
	_momentp = mp(1);
	_momentv = mv(1);
	_clp = _liftp/(qinf*sref);
	_clv = _liftv/(qinf*sref);
	_cdp = _dragp/(qinf*sref);
	_cdv = _dragv/(qinf*sref);
	_cmp = _momentp/(qinf*sref*lref);
	_cmv = _momentv/(qinf*sref*lref);
	
	// Compute section forces and moments
	
#pragma omp parallel for private(i)
	for ( i = 0; i < _nspan; i++ )
	{
		_sections[i].computeForceMoment(alpha, uinf, rhoinf, viscous);
	}
}

double Wing::lift () const { return _liftp + _liftv; }
const double & Wing::pressureLift () const { return _liftp; }
const double & Wing::viscousLift () const { return _liftv; }

double Wing::drag () const { return _dragp + _dragv; }
const double & Wing::pressureDrag () const { return _dragp; }
const double & Wing::viscousDrag () const { return _dragv; }

double Wing::pitchingMoment () const { return _momentp + _momentv; }
const double & Wing::pressurePitchingMoment () const { return _momentp; }
const double & Wing::viscousPitchingMoment () const { return _momentv; }

double Wing::liftCoefficient () const { return _clp + _clv; }
const double & Wing::pressureLiftCoefficient () const { return _clp; }
const double & Wing::viscousLiftCoefficient () const { return _clv; }

double Wing::dragCoefficient () const { return _cdp + _cdv; }
const double & Wing::pressureDragCoefficient () const { return _cdp; }
const double & Wing::viscousDragCoefficient () const { return _cdv; }

double Wing::pitchingMomentCoefficient () const { return _cmp + _cmv; }
const double & Wing::pressurePitchingMomentCoefficient () const { return _cmp; }
const double & Wing::viscousPitchingMomentCoefficient () const { return _cmv; }

/******************************************************************************/
//
// Writes forces and moments for aircraft and each wing to CSV formatted files
//
/******************************************************************************/
int Wing::writeForceMoment ( int iter ) const
{
  std::ofstream f;
  std::string fname;

  fname = "forcemoment/" + _name + "_forcemoment.csv";

  // Write header during first iteration

  if (iter == 1)
  {
    f.open(fname.c_str(), std::fstream::out);
    if (! f.is_open())
    {
      print_warning("Wing::writeForceMoment",
                    "Unable to open " + fname + " for writing.");
      return 1;
    }
	if (viscous)
	{
      f << "\"Iter\",\"Lift\",\"Liftp\",\"Liftv\","
		<<          "\"Drag\",\"Dragp\",\"Dragv\","
		<<"\"Pitching_moment\",\"Pitching_momentp\",\"Pitching_momentv\","
		<<"\"CL\",\"CLp\",\"CLv\","
		<<"\"CD\",\"CDp\",\"CDv\","
		<<"\"Cm\",\"Cmp\",\"Cmv\"" << std::endl;
	}
	else
	{
      f << "\"Iter\",\"Lift\",\"Drag\",\"Pitching moment\","
        << "\"CL\",\"CD\",\"Cm\"" << std::endl;
	}
  }
  else
  {
    f.open(fname.c_str(), std::fstream::app);
    if (! f.is_open())
    {
      print_warning("Wing::writeForceMoment",
                    "Unable to open " + fname + " for writing.");
      return 1;
    }
  }
  f << iter << ",";
  f.setf(std::ios_base::scientific);
  f << std::setprecision(7);
  if (viscous)
  {
    f << _liftp + _liftv << ",";
    f << _liftp << ",";
    f << _liftv << ",";
    f << _dragp + _dragv << ",";
    f << _dragp << ",";
    f << _dragv << ",";
    f << _momentp + _momentv << ",";
    f << _momentp << ",";
    f << _momentv << ",";
    f << _clp + _clv << ",";
    f << _clp << ",";
    f << _clv << ",";
    f << _cdp + _cdv << ",";
    f << _cdp << ",";
    f << _cdv << ",";
    f << _cmp + _cmv << ",";
    f << _cmp << ",";
    f << _cmv << std::endl;
  }
  else
  {
    f << _liftp << ",";
    f << _dragp << ",";
    f << _momentp << ",";
    f << _clp << ",";
    f << _cdp << ",";
    f << _cmp << std::endl;
  }

  f.close();

  return 0;
}

/******************************************************************************/
//
// Writes sectional forces to file
//
/******************************************************************************/
int Wing::writeSectionForceMoment ( int iter ) const
{
  std::ofstream f;
  std::string fname;
  std::vector<double> y_flat;
  int i;
  double dy, dz, ds;

  // Compute "flattened" spanwise locations

  y_flat.resize(_nspan);
  y_flat[0] = 0.;
  for ( i = 1; i < int(_nspan); i++ )
  {
    dy = _sections[i].y() - _sections[i-1].y();
    dz = _sections[i].zle() - _sections[i-1].zle();
    ds = std::sqrt(dy*dy + dz*dz);
    y_flat[i] = y_flat[i-1] + ds;
  }

  fname = "sectional/" + _name + "_sectional_iter" + int2string(iter) + ".csv";

  // Write header

  f.open(fname.c_str(), std::fstream::out);
  if (! f.is_open())
  {
    print_warning("Wing::writeSectionForces",
                  "Unable to open " + fname + " for writing.");
    return 1;
  }
  f << "\"xle\",\"y\",\"y_flat\",\"zle\",\"c\",";
  if (viscous)
  {
	f << "\"Re\",\"Cl\",\"Clp\",\"Clv\",\"Cd\",\"Cdp\",\"Cdv\","
	  <<        "\"Cm\",\"Cmp\",\"Cmv\",\"cCl\",\"cCd\"" << std::endl;
  }
  else
  {
	f << "\"Cl\",\"Cd\",\"Cm\",\"cCl\",\"cCd\"" << std::endl;
  }

  // Write data for sections and mirror image
  
  f.setf(std::ios_base::scientific);
  f << std::setprecision(7);
  for ( i = _nspan-1; i >= 0; i-- )
  {
    f << _sections[i].xle() << ",";
    f << _sections[i].y() << ",";
    f << y_flat[i] << ",";
    f << _sections[i].zle() << ",";
    f << _sections[i].chord() << ",";
	if (viscous)
	{
	  f << _sections[i].reynoldsNumber() << ",";
	  f << _sections[i].liftCoefficient() << ",";
	  f << _sections[i].pressureLiftCoefficient() << ",";
	  f << _sections[i].viscousLiftCoefficient() << ",";
	  f << _sections[i].dragCoefficient() << ",";
	  f << _sections[i].pressureDragCoefficient() << ",";
	  f << _sections[i].viscousDragCoefficient() << ",";
	  f << _sections[i].pitchingMomentCoefficient() << ",";
	  f << _sections[i].pressurePitchingMomentCoefficient() << ",";
	  f << _sections[i].viscousPitchingMomentCoefficient() << ",";
	  f << _sections[i].chord() * _sections[i].liftCoefficient() << ",";
	  f << _sections[i].chord() * _sections[i].dragCoefficient() << std::endl;
	}
	else
	{
	  f << _sections[i].liftCoefficient() << ",";
	  f << _sections[i].dragCoefficient() << ",";
	  f << _sections[i].pitchingMomentCoefficient() << ",";
	  f << _sections[i].chord() * _sections[i].liftCoefficient() << ",";
	  f << _sections[i].chord() * _sections[i].dragCoefficient() << std::endl;
	}
  }
  for ( i = 1; i < int(_nspan); i++ )
  {
    f << _sections[i].xle() << ",";
    f << -_sections[i].y() << ",";
    f << -y_flat[i] << ",";
    f << _sections[i].zle() << ",";
    f << _sections[i].chord() << ",";
	if (viscous)
	{
	  f << _sections[i].reynoldsNumber() << ",";
	  f << _sections[i].liftCoefficient() << ",";
	  f << _sections[i].pressureLiftCoefficient() << ",";
	  f << _sections[i].viscousLiftCoefficient() << ",";
	  f << _sections[i].dragCoefficient() << ",";
	  f << _sections[i].pressureDragCoefficient() << ",";
	  f << _sections[i].viscousDragCoefficient() << ",";
	  f << _sections[i].pitchingMomentCoefficient() << ",";
	  f << _sections[i].pressurePitchingMomentCoefficient() << ",";
	  f << _sections[i].viscousPitchingMomentCoefficient() << ",";
	  f << _sections[i].chord() * _sections[i].liftCoefficient() << ",";
	  f << _sections[i].chord() * _sections[i].dragCoefficient() << std::endl;
	}
	else
	{
	  f << _sections[i].liftCoefficient() << ",";
	  f << _sections[i].dragCoefficient() << ",";
	  f << _sections[i].pitchingMomentCoefficient() << ",";
	  f << _sections[i].chord() * _sections[i].liftCoefficient() << ",";
	  f << _sections[i].chord() * _sections[i].dragCoefficient() << std::endl;
	}
  }
  f.close();

  return 0;
}
