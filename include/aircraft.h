// Header for aircraft class

#ifndef AIRCRAFT_H
#define AIRCRAFT_H

#include <vector>
#include <string>
#include "wing.h"

class Vertex;
class Face;

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

    // Set up pointers to vertices, panels, and wake elements

    void setGeometryPointers ();

  public:

    // Constructor

    Aircraft ();

    // Read from XML

    int readXML ( const std::string & geom_file );

    // Write VTK viz

    int writeViz ( const std::string & fname ) const;
};

#endif
