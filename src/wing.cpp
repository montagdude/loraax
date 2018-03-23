#define _USE_MATH_DEFINES

#include <vector>
#include <cmath>
#include <Eigen/Core>
#include "algorithms.h"
#include "util.h"
#include "sectional_object.h"
#include "airfoil.h"
#include "section.h"
#include "wing.h"

#include <iostream>

/******************************************************************************/
//
// Wing class. Contains sections, airfoils, paneling info, faces, a wake, etc.
// 
/******************************************************************************/

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

  // Sort airfoils

  nfoils = foils.size();
  if (nfoils == 1)
  {
    _foils = foils;
    return;
  }

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
  unsigned int i, j, nsecs;
  std::vector<double> sw, nom_stations;
  Eigen::Vector2d secvec;
  double a4, a5, rootsp, tipsp;

  if (_foils.size() < 1)
  {
    conditional_stop(1, "Wing::setupSections",
                     "At least one airfoil is required.");
    return 1;
  }

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
    return 2;
  }
#endif

  // Compute wing length and roll angles

  sw.resize(nsecs);
  sw[0] = 0.;
  for ( i = 1; i < nsecs; i++ )
  {
    secvec[0] = user_sections[i].y() - user_sections[i-1].y();
    secvec[1] = user_sections[i].zle() - user_sections[i-1].zle();
    sw[i] = sw[i-1] + secvec.norm();
    if (sw[i] == sw[i-1])
    {
      conditional_stop(1, "Wing::setupSections",
                       "Wing has duplicate Y sections.");
      return 2;
    }
    user_sections[i].setRoll(tan(secvec[1]/secvec[0])*180./M_PI);
  }

  // Compute nominal discretized section locations (stations) from spacing

  rootsp = _rootsprat * sw[nsecs-1] / double(_nspan-1);
  tipsp = _tipsprat * sw[nsecs-1] / double(_nspan-1);
  opt_tanh_spacing(_nspan, sw[nsecs-1], rootsp, tipsp, a4, a5);
  nom_stations.resize(_nspan);
  nom_stations[0] = 0.; 
  for ( i = 1; i < _nspan; i++ )
  {
    nom_stations[i] = nom_stations[i-1] +
      tanh_spacing(i-1, a4, a5, _nspan, sw[nsecs-1], rootsp, tipsp);
  }

  // Optimize stations to be as close as possible to nominal but lining up with
  // user sections

  return 0; 
}
