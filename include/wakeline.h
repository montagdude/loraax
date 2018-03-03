// Header for WakeLine class

#ifndef WAKELINE_H
#define WAKELINE_H

#include <vector>
#include <Eigen/Core>

// Forward declarations

class Node;

/******************************************************************************/
//
// Wake line class.  Represents a collection of wake filament trailers behind
// a trailing edge node.
//
/******************************************************************************/
class WakeLine {

  private:

    unsigned int _numfilaments;  // Number of vortex filaments

    std::vector<double> _x, _y, _z; 
                                 // Coordinates of filament endpoints
    std::vector<double> _xnew, _ynew, _znew;
                                 // New coordinates used for relaxing the wake

    Node * _tenode;              // Trailing edge node that it is shed from

  public:

    // Constructor

    WakeLine ();

    // Setting and accessing number of filaments

    void setNumFilaments ( unsigned int );
    unsigned int numFilaments () const;

    // Setting and accessing filament endpoint coordinates

    void setEndpoint ( unsigned int, const double &, const double &, 
                       const double & );
    const double & x ( unsigned int ) const;
    const double & y ( unsigned int ) const;
    const double & z ( unsigned int ) const;

    // Setting new coordinates for wake relaxation

    void setNewEndpoint ( unsigned int, const double &, const double &,
                          const double & );
    void updateEndpoints ();

    // Setting and accessing TE node

    void setTENode ( Node * );
    Node & TENode () const;

    // Induced velocity coefficient at a point

    Eigen::Vector3d VCoeff ( const double &, const double &, const double &,
                             const double & ) const;
};

#endif
