#define _USE_MATH_DEFINES

#include <vector>
#include <cmath>
#include "util.h"
#include "settings.h"
#include "vertex.h"
#include "panel.h"
#include "tripanel.h"
#include "quadpanel.h"
#include "wake.h"

/******************************************************************************/
//
// Wake class.
//
/******************************************************************************/

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
Wake::Wake ()
{
	_nspan = 0;
	_nstream = 0;
	_verts.resize(0);
	_newx.resize(0);
	_newy.resize(0);
	_newz.resize(0);
	_topteverts.resize(0);
	_botteverts.resize(0);
	_tris.resize(0); 
	_quads.resize(0); 
}

/******************************************************************************/
//
// initialize with a vector of vertices along the wing trailing edge
//
/******************************************************************************/
void Wake::initialize ( const std::vector<Vertex *> & topteverts,
                        const std::vector<Vertex *> & botteverts,
                        int & next_global_vertidx, int & next_global_elemidx )
{
	int ntris, nquads;
	unsigned int i, j, blvert, brvert, trvert, tlvert;
	double x, y, z, xviz, yviz, zviz;
	
	_topteverts = topteverts;
	_botteverts = botteverts;
	
	// Determine number of rows to store based on inputs
	
	_nspan = _topteverts.size();
	_nstream = wakesteps+1;		// Number of relaxable streamwise 
								// points. One more will be added for
								// the trailing "infinite" panel.
	_verts.resize(_nspan*(_nstream+1));
	_newx.resize(_nspan*(_nstream-1));
	_newy.resize(_nspan*(_nstream-1));
	_newz.resize(_nspan*(_nstream-1));

	// Add vertices along freestream direction
	
	for ( i = 0; int(i) < _nspan; i++ )
	{
		for ( j = 0; int(j) < _nstream+1; j++ )
		{
			if (j == 0)
			{
				x = _topteverts[i]->x();
				y = _topteverts[i]->y();
				z = _topteverts[i]->z();
				_verts[i*(_nstream+1)].setCoordinates(x, y, z);
			}
			else if (int(j) < _nstream)
			{
				// Vertices that will roll up

				x = _verts[i*(_nstream+1)].x()
				  + std::cos(wakeangle*M_PI/180.)*uinf*double(j)*dt;
				y = _verts[i*(_nstream+1)].y();
				z = _verts[i*(_nstream+1)].z()
				  + std::sin(wakeangle*M_PI/180.)*uinf*double(j)*dt;
				_verts[i*(_nstream+1)+j].setCoordinates(x, y, z);
			}
			else
			{
				// Trailing vertex at "infinity"

				x = _verts[i*(_nstream+1)+_nstream-1].x()
				  + uinfvec(0)*1000.*rollupdist/uinf;
				y = _verts[i*(_nstream+1)+_nstream-1].y();
				z = _verts[i*(_nstream+1)+_nstream-1].z()
				  + uinfvec(2)*1000.*rollupdist/uinf;
				_verts[i*(_nstream+1)+j].setCoordinates(x, y, z);
				xviz = _verts[i*(_nstream+1)+_nstream-1].x() + uinfvec(0)*dt;
				yviz = _verts[i*(_nstream+1)+_nstream-1].y();
				zviz = _verts[i*(_nstream+1)+_nstream-1].z() + uinfvec(2)*dt;
				_verts[i*(_nstream+1)+j].setVizCoordinates(xviz, yviz, zviz);
			}
			_verts[i*(_nstream+1)+j].setIdx(next_global_vertidx);
			_verts[i*(_nstream+1)+j].setWakeTime(-double(j)*dt);
			next_global_vertidx += 1;
		}
	}

	// Create tri panels (which can roll up). Tri panels are used for this
	// portion of the wake because panels are required to be planar, and
	// with quads there would be significant distortion.

	ntris = (_nspan-1)*(_nstream-1)*2;
	_tris.resize(ntris);
	for ( i = 0; int(i) < _nspan-1; i++ )
	{
		for ( j = 0; int(j) < _nstream-1; j++ )
		{
			blvert = i*(_nstream+1)+j+1;
			brvert = (i+1)*(_nstream+1)+j+1;
			trvert = (i+1)*(_nstream+1)+j;
			tlvert = i*(_nstream+1)+j;
			
			_tris[i*(_nstream-1)*2+j*2].setIdx(next_global_elemidx);
			_tris[i*(_nstream-1)*2+j*2].addVertex(&_verts[tlvert]);
			_tris[i*(_nstream-1)*2+j*2].addVertex(&_verts[blvert]);
			_tris[i*(_nstream-1)*2+j*2].addVertex(&_verts[brvert]);
			next_global_elemidx += 1;
			
			_tris[i*(_nstream-1)*2+j*2+1].setIdx(next_global_elemidx);
			_tris[i*(_nstream-1)*2+j*2+1].addVertex(&_verts[brvert]);
			_tris[i*(_nstream-1)*2+j*2+1].addVertex(&_verts[trvert]);
			_tris[i*(_nstream-1)*2+j*2+1].addVertex(&_verts[tlvert]);
			next_global_elemidx += 1;
		}
	}

	// Create quad panels which extend to 1000*rollupdist
	
	nquads = _nspan-1;
	_quads.resize(nquads);
	for ( i = 0; int(i) < _nspan-1; i++ )
	{
		blvert = i*(_nstream+1)+_nstream;
		tlvert = i*(_nstream+1)+_nstream-1;
		trvert = (i+1)*(_nstream+1)+_nstream-1;
		brvert = (i+1)*(_nstream+1)+_nstream;
		_quads[i].setIdx(next_global_elemidx);
		_quads[i].addVertex(&_verts[tlvert]);
		_quads[i].addVertex(&_verts[blvert]);
		_quads[i].addVertex(&_verts[brvert]);
		_quads[i].addVertex(&_verts[trvert]);
		next_global_elemidx += 1;
	}
}

/******************************************************************************/
//
// Convects wake vertices downstream (a.k.a. wake rollup) using 1st-order
// Euler integration.
//
// Note that currently, the wake rollup procedure neglects any influence of
// wake panel source strengths, which are used to model the boundary layer in
// the wake for viscous cases. This is because the source influence becomes
// singular at the panel corners. This issue could be avoided with careful
// treatment of the "corner cases" (ha), but it is also probably not worth the
// effort.
//
// 1st-order integration is still pretty accurate here because the wake
// convection algorithm allows each point to take only 1 time step. Increasing
// the itegration order causes difficulty because the evaluation points end up
// very close to panel edges, but not on them (which would make the velocity
// contribution identically 0 due to the vortex core model). The integration
// order can potentially be increased if that problem can be worked around, but
// it's probably not worthwhile since this method is already pretty accurate.
//
/******************************************************************************/
void Wake::convectVertices ( const double & dt,
                             const std::vector<Panel *> & allsurf,
                             const std::vector<Panel *> & allwake )
{
	int i, j;
	unsigned int v, k, nsurfpan, nwakepan, nmove;
	Eigen::Vector3d k1, dvel;
	double x, y, z;
// FIXME: make this a function of smallest panel size
	const double rcore = 1.E-8;

	nsurfpan = allsurf.size();
	nwakepan = allwake.size();
	nmove = _nspan*(_nstream-1);

#pragma omp parallel for private(v,i,j,x,y,z,k,k1,dvel)
	for ( v = 0; v < nmove; v++ )
	{
		i = int(floor(v/(_nstream-1)));
		j = v - i*(_nstream-1);
		
		x = _verts[i*(_nstream+1)+j].x();
		y = _verts[i*(_nstream+1)+j].y();
		z = _verts[i*(_nstream+1)+j].z();
	
	  // Calculate velocity and update position
	
		if (j == 0)
		{
		  // At trailing edge, use average top + bottom surface velocity
		
			k1(0) = 0.5*(_topteverts[i]->data(2) + _botteverts[i]->data(2));
			k1(1) = 0.5*(_topteverts[i]->data(3) + _botteverts[i]->data(3));
			k1(2) = 0.5*(_topteverts[i]->data(4) + _botteverts[i]->data(4));
		}
		else
		{
		  // Elsewhere in the wake, sum the surface and wake influences
		
			k1 = uinfvec;
			for ( k = 0; k < nsurfpan; k++ )
			{
				dvel = allsurf[k]->inducedVelocity(x, y, z, false, "top", true);
				k1 += dvel;
			}
			for ( k = 0; k < nwakepan; k++ )
			{
				dvel = allwake[k]->vortexVelocity(x, y, z, rcore, true);
				k1 += dvel;
			}
		}
		if (i == 0)
			k1(1) = 0.;
		k1 *= dt;
		
		_newx[v] = x + k1(0);
		_newy[v] = y + k1(1);
		_newz[v] = z + k1(2);
	}
}

/******************************************************************************/
//
// Updates wake with new positions computed during the convect routine
//
/******************************************************************************/
void Wake::update ()
{
	unsigned int i, j, ntris, nquads;
	double x, y, z, xviz, yviz, zviz;
	
	// Update vertex positions. New position for vertex (i,j) is equal to
	// convected position of vertex (i,j-1).

#pragma omp parallel for private(i,j,x,y,z,xviz,yviz,zviz)
	for ( i = 0; int(i) < _nspan; i++ )
	{
		_verts[i*(_nstream+1)].incrementWakeTime(dt);
		for ( j = 1; int(j) < _nstream; j++ )
		{
			x = _newx[i*(_nstream-1)+j-1];
			y = _newy[i*(_nstream-1)+j-1];
			z = _newz[i*(_nstream-1)+j-1];
			_verts[i*(_nstream+1)+j].setCoordinates(x, y, z);
			_verts[i*(_nstream+1)+j].incrementWakeTime(dt);
		}
	
		// Trailing vertices at "infinity" extend along freestream direction
		
		x = _newx[i*(_nstream-1)+_nstream-2] + uinfvec(0)*1000.*rollupdist/uinf;
		y = _newy[i*(_nstream-1)+_nstream-2];
		z = _newz[i*(_nstream-1)+_nstream-2] + uinfvec(2)*1000.*rollupdist/uinf;
		_verts[i*(_nstream+1)+_nstream].setCoordinates(x, y, z);
		_verts[i*(_nstream+1)+_nstream].incrementWakeTime(dt);
		xviz = _newx[i*(_nstream-1)+_nstream-2] + uinfvec(0)*dt;
		yviz = _newy[i*(_nstream-1)+_nstream-2];
		zviz = _newz[i*(_nstream-1)+_nstream-2] + uinfvec(2)*dt;
		_verts[i*(_nstream+1)+_nstream].setVizCoordinates(xviz, yviz, zviz);
	}

	// Recompute panel geometry
	
	ntris = _tris.size(); 
#pragma omp parallel for private(i)
	for ( i = 0; i < ntris; i++ )
	{
		_tris[i].recomputeGeometry();
	}

	nquads = _quads.size();
#pragma omp parallel for private(i)
	for ( i = 0; i < nquads; i++ )
	{
		_quads[i].recomputeGeometry();
	}

	// Move doublets downstream

#pragma omp parallel for private(i,j)
	for ( i = 0; int(i) < _nspan-1; i++ )
	{
		_quads[i].setDoubletStrength(
		              _tris[i*(_nstream-1)*2+(_nstream-2)*2].doubletStrength());
		for ( j = _nstream-2; j >= 1; j-- )
		{
			_tris[i*(_nstream-1)*2+j*2].setDoubletStrength(
			                 _tris[i*(_nstream-1)*2+(j-1)*2].doubletStrength());
			_tris[i*(_nstream-1)*2+j*2+1].setDoubletStrength(
			               _tris[i*(_nstream-1)*2+(j-1)*2+1].doubletStrength());
		}
	}
}

/******************************************************************************/
//
// Interpolates mass defect from 2D BL solution to line of vertices (viscous
// cases only)
//
/******************************************************************************/
void Wake::interpMassDefectLine ( unsigned int idx,
                                  const std::vector<double> & sw,
                                  const std::vector<double> & dstarw,
                                  const std::vector<double> & uedgew )
{
	unsigned int j, k, nw, curr;
	double s, weight1, weight2; 

	nw = sw.size();
	curr = 0;
	s = 0.;
	for ( j = 0; int(j) < _nstream; j++ )
	{
		for ( k = curr; k < nw; k++ )
		{
			if ( (s >= sw[k]) && (s <= sw[k+1]) )
			{
				// Interpolate between BL wake vertices

				weight2 = (s - sw[k]) / (sw[k+1] - sw[k]);
				weight1 = 1. - weight2;
				_verts[idx*(_nstream+1)+j].setData(8,
					weight1*dstarw[k] + weight2*dstarw[k+1]);
				_verts[idx*(_nstream+1)+j].setData(10,
					weight1*uedgew[k] + weight2*uedgew[k+1]);
				curr = k;
				break;
			}
			else if (s > sw[nw-1])
			{
				// Beyond BL wake, leave deltastar and uedge constant. This will
				// result in source strength = 0 beyond the BL wake, which is
				// the appropriate limiting behavior.

				_verts[idx*(_nstream+1)+j].setData(8,
					_verts[idx*(_nstream+1)+j-1].data(8));
				_verts[idx*(_nstream+1)+j].setData(10,
					_verts[idx*(_nstream+1)+j-1].data(10));
				curr = nw-2;
				break;
			}
#ifdef DEBUG
			// This shouldn't happen, but put it here for debugging
			
			else if (s < sw[k])
				conditional_stop(1, "Wake::interpMassDefectLine",
				                 "Interpolation failed 1.");
#endif
		}
		s += _verts[idx*(_nstream+1)+j].distance(_verts[idx*(_nstream+1)+j+1]);
	}
	
	// Far downstream, displacement thickness is 0

	_verts[idx*(_nstream+1)+_nstream].setData(8, 0.0);
	_verts[idx*(_nstream+1)+_nstream].setData(10, uinf);
}

/******************************************************************************/
//
// Computes source strengths from mass defect (viscous cases only)
//
/******************************************************************************/
void Wake::computeSourceStrengths ()
{
	unsigned int i, j;
	double mdefl1, mdefl2, dmdefl, dsl;
	double mdefr1, mdefr2, dmdefr, dsr;

#pragma omp parallel for private(i,j,mdefl1,mdefl2,dsl,dmdefl,mdefr1,mdefr2,\
                                 dsr,dmdefr)
	for ( i = 0; int(i) < _nspan-1; i++ )
	{
		for ( j = 0; int(j) < _nstream; j++ )
		{
			mdefl1 = _verts[i*(_nstream+1)+j].data(8)
			       * _verts[i*(_nstream+1)+j].data(10);
			mdefl2 = _verts[i*(_nstream+1)+j+1].data(8)
			       * _verts[i*(_nstream+1)+j+1].data(10);
			dsl = _verts[i*(_nstream+1)+j].distance(_verts[i*(_nstream+1)+j+1]);
			dmdefl = (mdefl2 - mdefl1) / dsl;

			mdefr1 = _verts[(i+1)*(_nstream+1)+j].data(8)
			       * _verts[(i+1)*(_nstream+1)+j].data(10);
			mdefr2 = _verts[(i+1)*(_nstream+1)+j+1].data(8)
			       * _verts[(i+1)*(_nstream+1)+j+1].data(10);
			dsr = _verts[(i+1)*(_nstream+1)+j].distance(
			                    _verts[(i+1)*(_nstream+1)+j+1]);
			dmdefr = (mdefr2 - mdefr1) / dsr;

			if (int(j) < _nstream-1)
			{
				_tris[i*(_nstream-1)*2+j*2].setSourceStrength(
				                            0.67*dmdefl + 0.33*dmdefr);
				_tris[i*(_nstream-1)*2+j*2+1].setSourceStrength(
				                            0.33*dmdefl + 0.67*dmdefr);
			}
			else
				_quads[i].setSourceStrength(0.5*(dmdefl + dmdefr));
		}
	}
}

/******************************************************************************/
//
// Access to verts and panels
//
/******************************************************************************/
unsigned int Wake::nVerts () const { return _verts.size(); }
unsigned int Wake::nTris () const { return _tris.size(); }
unsigned int Wake::nQuads () const { return _quads.size(); }
Vertex * Wake::vert ( unsigned int vidx )
{
#ifdef DEBUG
  if (vidx >= _verts.size())
    conditional_stop(1, "Wake::vert", "Index out of range.");
#endif

  return &_verts[vidx];
}

TriPanel * Wake::triPanel ( unsigned int tidx )
{
#ifdef DEBUG
  if (tidx >= _tris.size())
    conditional_stop(1, "Wake::triPanel", "Index out of range.");
#endif

  return &_tris[tidx];
}

QuadPanel * Wake::quadPanel ( unsigned int qidx )
{
#ifdef DEBUG
  if (qidx >= _quads.size())
    conditional_stop(1, "Wake::quadPanel", "Index out of range.");
#endif

  return &_quads[qidx];
}
