#include <vector>
#include "util.h"
#include "airfoil.h"
#include "section.h"
#include "wing.h"

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
void Wing::setAirfoils ( const std::vector<Airfoil> & foils )
{
  //FIXME: sort and do other checks
  _foils = foils;
}

/******************************************************************************/
//
// Set sections based on user section inputs and spacing. User sections define
// the planform shape, and the _sections variable defines the final
// discretization.
//
/******************************************************************************/
int Wing::setupSections ( const std::vector<Section> & user_sections )
{
  if (_foils.size() < 1)
  {
    conditional_stop(1, "Wing::setupSections",
                     "At least one airfoil is required.");
    return 1;
  }

  return 0; 
}
