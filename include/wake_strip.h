// Header for WakeStrip class

#ifndef WAKESTRIP_H
#define WAKESTRIP_H

#include <vector>

class Panel;

/******************************************************************************/
//
// WakeStrip class. Stores pointers to connect trailing edge panels and the
// wake panels that they shed into, for the purpose of setting BCs correctly
// in the linear system.
//
/******************************************************************************/
class WakeStrip {

	private:

	Panel * _toptepan, * _bottepan;	// TE panels
	std::vector<Panel *> _panels;	// Strip of trailing vortex elements

	public:

	// Constructor
	
	WakeStrip ();
	
	// Add panels in a strip
	
	void setNPanels ( unsigned int npanels );
	int setPanelPointer ( unsigned int pidx, Panel * panel );
	
	// Set trailing edge panel pointers
	
	void setTEPanels ( Panel * toptepan, Panel * bottepan );
	
	// Access panels
	
	unsigned int nPanels () const;
	Panel * topTEPan ();
	Panel * botTEPan ();
	Panel * panel ( unsigned int pidx );
};

#endif
