// Header for WakeStrip class

#ifndef WAKESTRIP_H
#define WAKESTRIP_H

#include <vector>

class Panel;
class Vortex;

/******************************************************************************/
//
// WakeStrip class. Stores pointers to connect trailing edge panels and the
// wake vortices that they shed into, for the purpose of setting BCs correctly
// in the linear system.
//
/******************************************************************************/
class WakeStrip {

  private:

    Panel * _toptepan, * _bottepan;	// TE panels
    std::vector<Vortex *> _vorts;	// Strip of trailing vortex elements

  public:

    // Constructor

    WakeStrip ();

    // Add vortex elements in a strip

    void setNVortices ( unsigned int nvorts );
    int setVortexPointer ( unsigned int vidx, Vortex * vort );

    // Set trailing edge panel pointers

    void setTEPanels ( Panel * toptepan, Panel * bottepan );

    // Access panels and vortices

    unsigned int nVortices () const;
    Panel * topTEPan ();
    Panel * botTEPan ();
    Vortex * vortex ( unsigned int vidx );
};

#endif
