// Header for aircraft class

#ifndef AIRCRAFT_H
#define AIRCRAFT_H

#include <vector>
#include <string>
#include "wing.h"

class Vertex;
class Panel;
class Vortex;

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
    std::vector<Panel *> _panels;	// Pointers to panels 
    std::vector<Vertex *> _wakeverts;   // Pointers to wake vertices
    std::vector<Vortex *> _vorts;	// Pointers to wake vortex panels and
                                        //   horseshoes

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
