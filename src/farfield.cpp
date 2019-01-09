#include <vector>
#include <cmath>
#include "vertex.h"
#include "quadpanel.h"
#include "farfield.h"

/******************************************************************************/
//
// Farfield class.
//
/******************************************************************************/

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
Farfield::Farfield ()
{
    _nx = 0;
    _ny = 0;
    _nz = 0;
    _lenx = 0.;
    _leny = 0.;
    _lenz = 0.;
    _vertarray.resize(0);
    _quadarray.resize(0);
    _verts.resize(0);
    _quads.resize(0);
}

/******************************************************************************/
//
// Create farfield box
//
/******************************************************************************/
void Farfield::initialize ( unsigned int nx, unsigned int ny, unsigned int nz,
                            const double & lenx, const double & leny,
                            const double & lenz, const double & minf,
                            int & next_global_vertidx,
                            int & next_global_elemidx )
{
    unsigned int i, j;
    double beta, x, y, z, dx, dy, dz;

    _nx = nx;
    _ny = ny;
    _nz = nz;
    dx = _lenx / double(_nx-1);
    dy = _leny / double(_ny-1);
    dz = _lenz / double(_nz-1);
    beta = std::sqrt(1. - std::pow(minf, 2.));

    // Set up arrays. Each array has indices face, i, j. There are 6 faces, and
    // i and j go from 0 to _nx, _ny, or _nz.

    _vertarray.resize(6);
    _quadarray.resize(6);

    // Left face

    _vertarray[0].resize(_nx);
    y = -_leny/2.;
    for ( i = 0; i < _nx; i++ )
    {
        x = -_lenx/2. + double(i)*dx;
        _vertarray[0][i].resize(_nz);
        for ( j = 0; j < _nz; j++ )
        {
            z = _lenz/2. - double(j)*dz;
            _vertarray[0][i][j].setIdx(next_global_vertidx);
            _vertarray[0][i][j].setCoordinates(x, y, z);
            _vertarray[0][i][j].setIncompressibleCoordinates(x/beta, y, z);
            _verts.push_back(&_vertarray[0][i][j]);
            next_global_vertidx += 1;
        }
    }
    _quadarray[0].resize(_nx-1);
    for ( i = 0; i < _nx-1; i++ )
    {
        _quadarray[0][i].resize(_nz-1);
        for ( j = 0; j < _nz-1; j++ )
        {
            _quadarray[0][i][j].setIdx(next_global_elemidx);
            _quadarray[0][i][j].addVertex(&_vertarray[0][i][j]);
            _quadarray[0][i][j].addVertex(&_vertarray[0][i+1][j]);
            _quadarray[0][i][j].addVertex(&_vertarray[0][i+1][j+1]);
            _quadarray[0][i][j].addVertex(&_vertarray[0][i][j+1]);
            _quads.push_back(&_quadarray[0][i][j]);
            next_global_elemidx += 1;
        }
    }

    // Right face

    _vertarray[1].resize(_nx);
    y = _leny/2.;
    for ( i = 0; i < _nx; i++ )
    {
        x = -_lenx/2. + double(i)*dx;
        _vertarray[1][i].resize(_nz);
        for ( j = 0; j < _nz; j++ )
        {
            z = _lenz/2. - double(j)*dz;
            _vertarray[1][i][j].setIdx(next_global_vertidx);
            _vertarray[1][i][j].setCoordinates(x, y, z);
            _vertarray[1][i][j].setIncompressibleCoordinates(x/beta, y, z);
            _verts.push_back(&_vertarray[1][i][j]);
            next_global_vertidx += 1;
        }
    }
    _quadarray[1].resize(_nx-1);
    for ( i = 0; i < _nx-1; i++ )
    {
        _quadarray[1][i].resize(_nz-1);
        for ( j = 0; j < _nz-1; j++ )
        {
            _quadarray[1][i][j].setIdx(next_global_elemidx);
            _quadarray[1][i][j].addVertex(&_vertarray[1][i][j]);
            _quadarray[1][i][j].addVertex(&_vertarray[1][i][j+1]);
            _quadarray[1][i][j].addVertex(&_vertarray[1][i+1][j+1]);
            _quadarray[1][i][j].addVertex(&_vertarray[1][i+1][j]);
            _quads.push_back(&_quadarray[1][i][j]);
            next_global_elemidx += 1;
        }
    }

    // Front face

    _vertarray[2].resize(_ny);
    x = -_lenx/2.;
    for ( i = 0; i < _ny; i++ )
    {
        y = -_leny/2. + double(i)*dy;
        _vertarray[2][i].resize(_nz);
        for ( j = 0; j < _nz; j++ )
        {
            z = _lenz/2. - double(j)*dz;
            _vertarray[2][i][j].setIdx(next_global_vertidx);
            _vertarray[2][i][j].setCoordinates(x, y, z);
            _vertarray[2][i][j].setIncompressibleCoordinates(x/beta, y, z);
            _verts.push_back(&_vertarray[2][i][j]);
            next_global_vertidx += 1;
        }
    }
    _quadarray[2].resize(_ny-1);
    for ( i = 0; i < _ny-1; i++ )
    {
        _quadarray[2][i].resize(_nz-1);
        for ( j = 0; j < _nz-1; j++ )
        {
            _quadarray[2][i][j].setIdx(next_global_elemidx);
            _quadarray[2][i][j].addVertex(&_vertarray[2][i][j]);
            _quadarray[2][i][j].addVertex(&_vertarray[2][i][j+1]);
            _quadarray[2][i][j].addVertex(&_vertarray[2][i+1][j+1]);
            _quadarray[2][i][j].addVertex(&_vertarray[2][i+1][j]);
            _quads.push_back(&_quadarray[2][i][j]);
            next_global_elemidx += 1;
        }
    }

    // Back face

    _vertarray[3].resize(_ny);
    x = _lenx/2.;
    for ( i = 0; i < _ny; i++ )
    {
        y = -_leny/2. + double(i)*dy;
        _vertarray[3][i].resize(_nz);
        for ( j = 0; j < _nz; j++ )
        {
            z = _lenz/2. - double(j)*dz;
            _vertarray[3][i][j].setIdx(next_global_vertidx);
            _vertarray[3][i][j].setCoordinates(x, y, z);
            _vertarray[3][i][j].setIncompressibleCoordinates(x/beta, y, z);
            _verts.push_back(&_vertarray[3][i][j]);
            next_global_vertidx += 1;
        }
    }
    _quadarray[3].resize(_ny-1);
    for ( i = 0; i < _ny-1; i++ )
    {
        _quadarray[3][i].resize(_nz-1);
        for ( j = 0; j < _nz-1; j++ )
        {
            _quadarray[3][i][j].setIdx(next_global_elemidx);
            _quadarray[3][i][j].addVertex(&_vertarray[3][i][j]);
            _quadarray[3][i][j].addVertex(&_vertarray[3][i+1][j]);
            _quadarray[3][i][j].addVertex(&_vertarray[3][i+1][j+1]);
            _quadarray[3][i][j].addVertex(&_vertarray[3][i][j+1]);
            _quads.push_back(&_quadarray[3][i][j]);
            next_global_elemidx += 1;
        }
    }

    // Top face

    _vertarray[4].resize(_nx);
    z = _lenz/2.;
    for ( i = 0; i < _nx; i++ )
    {
        x = -_leny/2. + double(i)*dx;
        _vertarray[4][i].resize(_ny);
        for ( j = 0; j < _ny; j++ )
        {
            y = -_leny/2. + double(j)*dy;
            _vertarray[4][i][j].setIdx(next_global_vertidx);
            _vertarray[4][i][j].setCoordinates(x, y, z);
            _vertarray[4][i][j].setIncompressibleCoordinates(x/beta, y, z);
            _verts.push_back(&_vertarray[4][i][j]);
            next_global_vertidx += 1;
        }
    }
    _quadarray[4].resize(_nx-1);
    for ( i = 0; i < _nx-1; i++ )
    {
        _quadarray[4][i].resize(_ny-1);
        for ( j = 0; j < _ny-1; j++ )
        {
            _quadarray[4][i][j].setIdx(next_global_elemidx);
            _quadarray[4][i][j].addVertex(&_vertarray[4][i][j]);
            _quadarray[4][i][j].addVertex(&_vertarray[4][i][j+1]);
            _quadarray[4][i][j].addVertex(&_vertarray[4][i+1][j+1]);
            _quadarray[4][i][j].addVertex(&_vertarray[4][i+1][j]);
            _quads.push_back(&_quadarray[4][i][j]);
            next_global_elemidx += 1;
        }
    }

    // Bottom face

    _vertarray[5].resize(_nx);
    z = -_lenz/2.;
    for ( i = 0; i < _nx; i++ )
    {
        x = -_leny/2. + double(i)*dx;
        _vertarray[5][i].resize(_ny);
        for ( j = 0; j < _ny; j++ )
        {
            y = -_leny/2. + double(j)*dy;
            _vertarray[5][i][j].setIdx(next_global_vertidx);
            _vertarray[5][i][j].setCoordinates(x, y, z);
            _vertarray[5][i][j].setIncompressibleCoordinates(x/beta, y, z);
            _verts.push_back(&_vertarray[5][i][j]);
            next_global_vertidx += 1;
        }
    }
    _quadarray[5].resize(_nx-1);
    for ( i = 0; i < _nx-1; i++ )
    {
        _quadarray[5][i].resize(_ny-1);
        for ( j = 0; j < _ny-1; j++ )
        {
            _quadarray[5][i][j].setIdx(next_global_elemidx);
            _quadarray[5][i][j].addVertex(&_vertarray[5][i][j]);
            _quadarray[5][i][j].addVertex(&_vertarray[5][i+1][j]);
            _quadarray[5][i][j].addVertex(&_vertarray[5][i+1][j+1]);
            _quadarray[5][i][j].addVertex(&_vertarray[5][i][j+1]);
            _quads.push_back(&_quadarray[5][i][j]);
            next_global_elemidx += 1;
        }
    }
}
