// Header for VortexRing class

#ifndef VORTEXRING_H
#define VORTEXRING_H

#include <vector>
#include <Eigen/Core>

class Vertex;

/******************************************************************************/
//
// Vortex ring class
//
/******************************************************************************/
class VortexRing {

  private:

    int _idx;				// Identifier index
    unsigned int _currverts;		// Number of vertices
    double _gamma;			// Circulation strength
    std::vector<Vertex *> _verts;	// 4 endpoint vertices

  public:

    // Constructor

    VortexRing ();

    // Set / access index

    void setIdx ( int idx );
    int idx () const;

    // Set / access vertices

    int addVertex ( Vertex * vert );
    Vertex & vertex ( unsigned int vidx ) const;

    // Setting and accessing circulation

    void setCirculation ( const double & gamma );
    const double & circulation () const;

    // Induced velocity coefficient and induced velocity at a point

    Eigen::Vector3d VCoeff ( const double & x, const double & y,
                             const double & z, const double & rcore ) const;

    Eigen::Vector3d inducedVelocity ( const double & x, const double & y,
                                      const double & z,
                                      const double & rcore ) const;
};

#endif
