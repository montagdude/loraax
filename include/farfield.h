// Header for Farfield class

#ifndef FARFIELD_H
#define FARFIELD_H

#include <vector>
#include "vertex.h"
#include "quadpanel.h"

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

    public:

    // Constructor
    
    Farfield ();
    
    // Initialize
    
    void initialize ( unsigned int nx, unsigned int ny, unsigned int nz,
                      const double & lenx, const double & leny,
                      const double & lenz, int & next_global_vertidx,
                      int & next_global_elemidx );

    // Access vertices and panels
    
    unsigned int nVerts () const;
    unsigned int nQuads () const;
    Vertex * vert ( unsigned int vidx );
    QuadPanel * quadPanel ( unsigned int qidx );

};

#endif
