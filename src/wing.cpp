#define _USE_MATH_DEFINES

#include <vector>
#include <cmath>
#include <Eigen/Core>
#include "algorithms.h"
#include "util.h"
#include "settings.h"
#include "sectional_object.h"
#include "airfoil.h"
#include "section.h"
#include "wing.h"

#include <iostream>
#include <iomanip>

// Data for optimizing spanwise spacing

std::vector<double> nom_stations, new_stations, s_wing;
std::vector<int> fixed_stations;

/******************************************************************************/
//
// Wing class. Contains sections, airfoils, paneling info, faces, a wake, etc.
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
    dist2 += pow(s_wing[i+1]-nom_stations[comb[i]], 2.);
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
      objval += 5.*pow((space-nom_space)/nom_space, 2.);
    else
      objval += pow((space-nom_space)/nom_space, 2.);
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
// Default constructor
//
/******************************************************************************/
Wing::Wing ()
{
  _name = "";
  _nchord = 0;
  _nspan = 0;
  _lesprat = 1.;
  _tesprat = 1.;
  _rootsprat = 1.;
  _tipsprat = 1.;
  _sections.resize(0);
  _foils.resize(0);
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
    sorted_user_sections[i].setRoll(tan(secvec[1]/secvec[0])*180./M_PI);
  }

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
    _sections[i].airfoil().smoothPaneling(xfoil_geom_opts);

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
        roll = sorted_user_sections[j+1].roll();
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
  }

  return 0; 
}
