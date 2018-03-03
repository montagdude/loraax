// Header for Wake class

#ifndef WAKE_H
#define WAKE_H

#include <vector>
#include <Eigen/Core>
#include "wakeline.h"

// Forward declarations

class Node;

/******************************************************************************/
//
// Wake class. Stores wake lines and horseshoe vortices and facilitates their
// setup and connections with each other and the trailing edges of the aircraft
//
/******************************************************************************/
class Wake {

  private:

    unsigned int _numlines;      // Number of wake lines
    unsigned int _numhorseshoes; // Number of horseshoe vortices

    std::vector<WakeLine> _wakelines;
                                 // Vector of wake lines

  public:

    // Constructor

    Wake ();

    // Returns number of wake lines

    unsigned int numWakeLines () const;

    // Creates a wake line with given length and number of filaments, aligned
    // with the freestream

    void addWakeLine ( const double &, unsigned int, const Eigen::Vector3d &,
                       Node * ); 

    // Returns a reference to a wake line

    WakeLine & wakeLine ( unsigned int );
};

#endif
