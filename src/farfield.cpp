#include <vector>
#include <cmath>
#include <Eigen/Core>
#include "util.h"
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
    unsigned int i, j, k;
    double beta, x, y, z, dx, dy, dz;

    _nx = nx;
    _ny = ny;
    _nz = nz;
    _lenx = lenx;
    _leny = leny;
    _lenz = lenz;
    dx = _lenx / double(_nx-1);
    dy = _leny / double(_ny-1);
    dz = _lenz / double(_nz-1);
    beta = std::sqrt(1. - std::pow(minf, 2.));

    // Set up arrays. Each array has indices face, i, j. There are 6 faces, and
    // i and j go from 0 to _nx, _ny, or _nz.

    _vertarray.resize(6);
    _quadarray.resize(6);

    // Left face. Note faces are set up with normals pointing out of the volume,
    // as needed for the fluid rate of change of momentum calculation.

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
            _quadarray[0][i][j].addVertex(&_vertarray[0][i][j+1]);
            _quadarray[0][i][j].addVertex(&_vertarray[0][i+1][j+1]);
            _quadarray[0][i][j].addVertex(&_vertarray[0][i+1][j]);
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
            _quadarray[1][i][j].addVertex(&_vertarray[1][i+1][j]);
            _quadarray[1][i][j].addVertex(&_vertarray[1][i+1][j+1]);
            _quadarray[1][i][j].addVertex(&_vertarray[1][i][j+1]);
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
            _quadarray[2][i][j].addVertex(&_vertarray[2][i+1][j]);
            _quadarray[2][i][j].addVertex(&_vertarray[2][i+1][j+1]);
            _quadarray[2][i][j].addVertex(&_vertarray[2][i][j+1]);
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
            _quadarray[3][i][j].addVertex(&_vertarray[3][i][j+1]);
            _quadarray[3][i][j].addVertex(&_vertarray[3][i+1][j+1]);
            _quadarray[3][i][j].addVertex(&_vertarray[3][i+1][j]);
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
            _quadarray[4][i][j].addVertex(&_vertarray[4][i+1][j]);
            _quadarray[4][i][j].addVertex(&_vertarray[4][i+1][j+1]);
            _quadarray[4][i][j].addVertex(&_vertarray[4][i][j+1]);
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
            _quadarray[5][i][j].addVertex(&_vertarray[5][i][j+1]);
            _quadarray[5][i][j].addVertex(&_vertarray[5][i+1][j+1]);
            _quadarray[5][i][j].addVertex(&_vertarray[5][i+1][j]);
            next_global_elemidx += 1;
        }
    }

    // Store pointers to vertices and quads

    for ( k = 0; k < 6; k++ )
    {
        for ( i = 0; i < _vertarray[k].size(); i++ )
        {
            for ( j = 0; j < _vertarray[k][i].size(); j++ )
            {
                _verts.push_back(&_vertarray[k][i][j]);
            }
        }
    }
    for ( k = 0; k < 6; k++ )
    {
        for ( i = 0; i < _quadarray[k].size(); i++ )
        {
            for ( j = 0; j < _quadarray[k][i].size(); j++ )
            {
                _quads.push_back(&_quadarray[k][i][j]);
            }
        }
    }
}

/*******************************************************************************

Access vertices and panels

*******************************************************************************/
unsigned int Farfield::nVerts () const { return _verts.size(); }
unsigned int Farfield::nQuads () const { return _quads.size(); }
Vertex * Farfield::vert ( unsigned int vidx )
{
#ifdef DEBUG
    if (vidx >= _verts.size())
        conditional_stop(1, "Farfield::vert", "Index out of range.");
#endif

    return _verts[vidx];
}

QuadPanel * Farfield::quadPanel ( unsigned int qidx )
{
#ifdef DEBUG
    if (qidx >= _quads.size())
        conditional_stop(1, "Farfield::quad", "Index out of range.");
#endif

    return _quads[qidx];
}

/******************************************************************************

Velocity and pressure computation. Induced velocities are computed in
incompressible coordinates. Compressible velocity is: V_c = (U_i/beta, V_i, W_i)

*******************************************************************************/
void Farfield::computeVelocity ( const Eigen::Vector3d & uinfvec,
                                 const double & minf,
                                 const std::vector<Panel *> & allsurf,
                                 const std::vector<Panel *> & allwake )
{
    unsigned int i, k, nverts, nsurfpan, nwakepan;
    double beta, xinc, yinc, zinc;
    Eigen::Vector3d vel, dvel, dvelcomp;

    beta = std::sqrt(1. - std::pow(minf, 2.));
    nsurfpan = allsurf.size();
    nwakepan = allwake.size();

    nverts = nVerts();
#pragma omp parallel for private(i,xinc,yinc,zinc,vel,k,dvel,dvelcomp)
    for ( i = 0; i < nverts; i++ )
    {
        xinc = _verts[i]->xInc();
        yinc = _verts[i]->yInc();
        zinc = _verts[i]->zInc();
        vel = uinfvec;
        for ( k = 0; k < nsurfpan; k++ )
        {
            dvel = allsurf[k]->inducedVelocity(xinc, yinc, zinc, false, "top",
                                               true);
            dvelcomp << dvel(0)/beta, dvel(1), dvel(2);
            vel += dvelcomp;
        }
        for ( k = 0; k < nwakepan; k++ )
        {
            dvel = allwake[k]->vortexVelocity(xinc, yinc, zinc, 1e-08, true);
            dvelcomp << dvel(0)/beta, dvel(1), dvel(2);
            vel += dvelcomp;
        }
        _verts[i]->setData(2, vel(0));
        _verts[i]->setData(3, vel(1));
        _verts[i]->setData(4, vel(2));
    }
}

int Farfield::computePressure ( const double & uinf, const double & rhoinf,
                                const double & pinf )
{
    int retval;
    unsigned int i, nverts;
    double uinf2, vel2, ainf2, qinf, cpinc, minf2, beta, p0, cp, p, m2, gamm1,
           gamma, mach, rho0, rho;
    Eigen::Vector3d vel;

    uinf2 = std::pow(uinf, 2.);
    qinf = 0.5*rhoinf*uinf2;
    ainf2 = 1.4*pinf/rhoinf;
    minf2 = uinf2/ainf2;
    beta = std::sqrt(1. - minf2);
    gamma = 1.4;
    gamm1 = gamma - 1.;
    retval = 0;

    nverts = nVerts();
#pragma omp parallel for private(i,vel,vel2,cpinc,cp,p,p0,m2,mach)
    for ( i = 0; i < nverts; i++ )
    {
        vel << _verts[i]->data(2), _verts[i]->data(3), _verts[i]->data(4);
        vel2 = vel.squaredNorm();
        cpinc = 1. - vel2 / uinf2;

        // Prandtl-Glauert compressibility correction for pressure
  
        cp = cpinc / beta;
        p = cp*qinf + pinf;
  
        // Physical limits on pressure
  
        p0 = pinf * std::pow(1. + 0.5*gamm1*minf2, gamma/gamm1);
        if (p > p0)
        {
            p = p0;
            cp = (p - pinf) / qinf;
        }
  
        // Use isentropic relations to get local Mach number
  
        m2 = 2./gamm1 * (std::pow(p0/p, gamm1/gamma) - 1.);
        if (m2 > 1.0)
        {
            conditional_stop(1, "Farfield::computePressure",
                             "Locally supersonic flow detected.");
            retval = 1;
        }
        mach = std::sqrt(m2);
        rho0 = rhoinf * std::pow(1. + 0.5*gamm1*minf2, 1./gamm1);
        rho = rho0 / std::pow(1. + 0.5*gamm1*m2, 1./gamm1);

        _verts[i]->setData(5, p);
        _verts[i]->setData(6, cp);
        _verts[i]->setData(7, mach);
        _verts[i]->setData(8, rho);
    }
  
    return retval;
}
