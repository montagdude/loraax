#include <vector>
#include "util.h"
#include "panel.h"
#include "wake_strip.h"

/******************************************************************************/
//
// WakeStrip class. Stores pointers to connect trailing edge panels and the wake
// panels that they shed into, for the purpose of setting BCs correctly in the
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
	_panels.resize(0);
}

/******************************************************************************/
//
// Add panels in a strip
//
/******************************************************************************/
void WakeStrip::setNPanels ( unsigned int npanels )
{
#ifdef DEBUG
	if (_panels.size() != 0)
		print_warning("WakeStrip::setNPanels", "Clearing existing wake strip.");
#endif

	_panels.clear();
	_panels.resize(npanels);
}

int WakeStrip::setPanelPointer ( unsigned int pidx, Panel * panel )
{
#ifdef DEBUG
	if (pidx >= _panels.size())
	{
		conditional_stop(1, "WakeStrip::setPanelPointer", "Index out of range.");
		return 1;
	}
#endif

	_panels[pidx] = panel;

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
unsigned int WakeStrip::nPanels () const { return _panels.size(); }
Panel * WakeStrip::topTEPan () { return _toptepan; }
Panel * WakeStrip::botTEPan () { return _bottepan; }
Panel * WakeStrip::panel ( unsigned int pidx )
{
#ifdef DEBUG
	if (pidx >= _panels.size())
		conditional_stop(1, "WakeStrip::panel", "Index out of range.");
#endif

	return _panels[pidx];
}
