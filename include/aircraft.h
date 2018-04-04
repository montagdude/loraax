// Header for aircraft class

#ifndef AIRCRAFT_H
#define AIRCRAFT_H

#include <vector>
#include <string>
#include "wing.h"

class Vertex;
class Face;
class VortexRing;

/******************************************************************************/
//
// Aircraft class. Contains some number of wings and related data and members.
//
/******************************************************************************/
class Aircraft {

  private:

    std::vector<Wing> _wings;		// Wings

    double _sref;			// Reference area
    double _lref;			// Pitching moment reference length
    Eigen::Vector3d _momcen;		// Pitching moment reference point

    std::vector<Vertex *> _verts;	// Pointers to vertices
    std::vector<Face *> _panels;	// Pointers to panels 
    std::vector<Vertex *> _wakeverts;   // Pointers to wake vertices
    std::vector<VortexRing *> _vrings;	// Pointers to vortex rings

    // Set up pointers to vertices, panels, and wake elements

    void setGeometryPointers ();

    // Write VTK viz

    int writeSurfaceViz ( const std::string & prefix ) const;
    int writeWakeViz ( const std::string & prefix ) const;

  public:

    // Constructor

    Aircraft ();

    // Read from XML

    int readXML ( const std::string & geom_file );

    // Write VTK viz

    int writeViz ( const std::string & prefix ) const;
};

#endif
