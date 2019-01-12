// Header for Farfield class

#ifndef FARFIELD_H
#define FARFIELD_H

#include <vector>
#include <Eigen/Core>
#include "vertex.h"
#include "quadpanel.h"

class Panel;

/******************************************************************************/
//
// Farfield class. Used for calculating rate of change of momentum of air in the
// farfield.
//
/******************************************************************************/
class Farfield {

    private:

    unsigned int _nx, _ny, _nz;         // Number of panels in each direction
    double _lenx, _leny, _lenz;         // Dimensions
    std::vector< std::vector< std::vector<Vertex> > > _vertarray;
    std::vector< std::vector< std::vector<QuadPanel> > > _quadarray;
    std::vector<Vertex *> _verts;
    std::vector<QuadPanel *> _quads;
    Eigen::Vector3d _momrate, _pforce, _aeroforce;
    double _lift, _induced_drag, _cl, _cdi;

    public:

    // Constructor
    
    Farfield ();
    
    // Initialize
    
    void initialize ( unsigned int nx, unsigned int ny, unsigned int nz,
                      const double & lenx, const double & leny,
                      const double & lenz, const double & minf,
                      int & next_global_vertidx,
                      int & next_global_elemidx );

    // Access vertices and panels
    
    unsigned int nVerts () const;
    unsigned int nQuads () const;
    Vertex * vert ( unsigned int vidx );
    QuadPanel * quadPanel ( unsigned int qidx );

    // Velocity and pressure (and cp, mach, and density) calculation

    void computeVelocity ( const Eigen::Vector3d & uinfvec, const double & minf,
                           const std::vector<Panel *> & allsurf,
                           const std::vector<Panel *> & allwake );
    int computePressure ( const double & uinf, const double & rhoinf,
                          const double & pinf );

    // Computes force (rate of change of momentum) on fluid inside control
    // volume. Also computes pressure force on fluid due to outer boundary.
    // Inviscid aero force acting on aircraft is pressureForce - force.

    void computeForce ( const double & alpha, const double & rhoinf,
                        const double & uinf, const double & sref );
    const Eigen::Vector3d & force () const;
    const Eigen::Vector3d & pressureForce () const;
};

#endif
