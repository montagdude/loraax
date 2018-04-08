#include <vector>
#include "util.h"
#include "panel.h"
#include "vortex.h"
#include "wake_strip.h"

/******************************************************************************/
//
// WakeStrip class. Stores pointers to connect trailing edge panels and the wake
// vortices that they shed into, for the purpose of setting BCs correctly in the
// linear system.
//
/******************************************************************************/

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
WakeStrip::WakeStrip ()
{
  _toptepan = NULL;
  _bottepan = NULL;
  _vorts.resize(0);
}

/******************************************************************************/
//
// Add vortex elements in a strip
//
/******************************************************************************/
void WakeStrip::setNVortices ( unsigned int nvorts )
{
#ifdef DEBUG
  if (_vorts.size() != 0)
    print_warning("WakeStrip::setNVortices", "Clearing existing wake strip.");
#endif

  _vorts.clear();
  _vorts.resize(nvorts);
}

int WakeStrip::setVortexPointer ( unsigned int vidx, Vortex * vort )
{
#ifdef DEBUG
  if (vidx >= _vorts.size())
  {
    conditional_stop(1, "WakeStrip::setVortexPointer", "Index out of range.");
    return 1;
  }
#endif

  _vorts[vidx] = vort;

  return 0;
} 

/******************************************************************************/
//
// Set trailing edge panel pointers
//
/******************************************************************************/
void WakeStrip::setTEPanels ( Panel * toptepan, Panel * bottepan )
{
#ifdef DEBUG
  if (_toptepan != NULL)
    print_warning("WakeStrip::setTEPanels",
                  "Overwriting existing top TE panel pointer.");

  if (_bottepan != NULL)
    print_warning("WakeStrip::setTEPanels",
                  "Overwriting existing bot TE panel pointer.");
#endif

  _toptepan = toptepan;
  _bottepan = bottepan;
}

/******************************************************************************/
//
// Access to panels and vortces
//
/******************************************************************************/
unsigned int WakeStrip::nVortices () const { return _vorts.size(); }
Panel * WakeStrip::topTEPan () { return _toptepan; }
Panel * WakeStrip::botTEPan () { return _bottepan; }
Vortex * WakeStrip::vortex ( unsigned int vidx )
{
#ifdef DEBUG
  if (vidx >= _vorts.size())
    conditional_stop(1, "WakeStrip::vortex", "Index out of range.");
#endif

  return _vorts[vidx];
}
