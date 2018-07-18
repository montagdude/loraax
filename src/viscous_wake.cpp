#include <vector>
#include "util.h"
#include "tripanel.h"
#include "vertex.h"
#include "section.h"
#include "viscous_wake.h"

/******************************************************************************/
//
// Viscous wake class. This is the non-convecting wake which provides the proper
// smooth decrease in displacement thickness behind the trailing edge, which is
// needed for correct pressure drag and smooth pressure distribution near the
// trailing edge in viscous cases. The viscous wake points are set by the
// boundary layer solution.
//
/******************************************************************************/

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
ViscousWake::ViscousWake ()
{
	_nspan = 0;
	_nwake = 0;
	_verts.resize(0);
	_tris.resize(0);
}

/******************************************************************************/
//
// Initializes viscous wake with sections from the wing. Note that the 2D
// section wake is not initialized until after the BL solution is computed for
// the first time, so this must be called after that point.
//
/******************************************************************************/
void ViscousWake::initialize ( std::vector<Section> & sections,
                               int & next_global_vertidx,
                               int & next_global_elemidx )
{
	unsigned int i, j, vcounter, tcounter;
	unsigned int blvert, brvert, trvert, tlvert;

	_nspan = sections.size();
	_nwake = sections[0].nWake();

	// Set up vertices

	_verts.resize(_nspan*_nwake);
	vcounter = 0;
	for ( i = 0; i < _nspan; i++ )
	{
#ifdef DEBUG
		if (sections[i].nWake() != _nwake)
		{
			conditional_stop(1, "ViscousWake::initialize",
			       "Sections must have the same number of points in the wake.");
		}
#endif
		for ( j = 0; j < _nwake; j++ )
		{
			_verts[vcounter] = &sections[i].wakeVert(j);
			_verts[vcounter]->setIdx(next_global_vertidx);
			next_global_vertidx += 1;
			vcounter += 1;
		}
	}

	// Set up tri panels

	_tris.resize((_nspan-1)*(_nwake-1)*2);
	tcounter = 0;
	for ( i = 0; i < _nspan-1; i++ )
	{
		for ( j = 0; j < _nwake-1; j++ )
		{
			blvert = i*_nwake+j+1;
			brvert = (i+1)*_nwake+j+1;
			trvert = (i+1)*_nwake+j;
			tlvert = i*_nwake+j;

			_tris[tcounter].setIdx(next_global_elemidx);
			_tris[tcounter].addVertex(_verts[tlvert]);
			_tris[tcounter].addVertex(_verts[blvert]);
			_tris[tcounter].addVertex(_verts[brvert]);
			next_global_elemidx += 1;
			tcounter += 1;

			_tris[tcounter].setIdx(next_global_elemidx);
			_tris[tcounter].addVertex(_verts[brvert]);
			_tris[tcounter].addVertex(_verts[trvert]);
			_tris[tcounter].addVertex(_verts[tlvert]);
			next_global_elemidx += 1;
			tcounter += 1;
		}
	}
}

/******************************************************************************/
//
// Updates panel geometry (vertices could have moved) and computes source
// strength from mass defect derivative
//
/******************************************************************************/
void ViscousWake::update ()
{
	unsigned int i, j, ntri;
	double dsl, mdefl1, mdefl2, dmdefl, dsr, mdefr1, mdefr2, dmdefr;
	Vertex *tlvert, *blvert, *trvert, *brvert;

	ntri = _tris.size();
#pragma omp parallel for private(i)
	for ( i = 0; i < ntri; i++ )
	{
		_tris[i].recomputeGeometry();
	}

	// Compute mass defect derivative and source strength

#pragma omp parallel for private(i,j,dsl,mdefl1,mdefl2,dmdefl,dsr,mdefr1,\
	                             mdefr2,dmdefr,tlvert,blvert,trvert,brvert)
	for ( i = 0; i < _nspan-1; i++ )
	{
		for ( j = 0; j < _nwake-1; j++ )
		{
			tlvert = _verts[i*_nwake+j];
			blvert = _verts[i*_nwake+j+1];
			dsl = tlvert->distance(*blvert);
			mdefl1 = tlvert->data(9)*tlvert->data(11);
			mdefl2 = blvert->data(9)*blvert->data(11);
			dmdefl = (mdefl2 - mdefl1) / dsl;

			trvert = _verts[(i+1)*_nwake+j];
			brvert = _verts[(i+1)*_nwake+j+1];
			dsr = trvert->distance(*brvert);
			mdefr1 = trvert->data(9)*trvert->data(11);
			mdefr2 = brvert->data(9)*brvert->data(11);
			dmdefr = (mdefr2 - mdefr1) / dsr;

			_tris[i*(_nwake-1)*2+j*2].setSourceStrength(0.67*dmdefl +
			                                            0.33*dmdefr);
			_tris[i*(_nwake-1)*2+j*2+1].setSourceStrength(0.33*dmdefl +
			                                              0.67*dmdefr);
		}
	}
}

/******************************************************************************/
//
// Access vertices and panels
//
/******************************************************************************/
unsigned int ViscousWake::nVerts () const { return _verts.size(); }
unsigned int ViscousWake::nTris () const { return _tris.size(); }
Vertex * ViscousWake::vert ( unsigned int vidx )
{
#ifdef DEBUG
	if (vidx >= _verts.size())
		conditional_stop(1, "ViscousWake::vert", "Index out of range.");
#endif

	return _verts[vidx];
}

TriPanel * ViscousWake::triPanel ( unsigned int tidx )
{
#ifdef DEBUG
	if (tidx >= _tris.size())
		conditional_stop(1, "ViscousWake::triPanel", "Index out of range.");
#endif

	return &_tris[tidx];
}
