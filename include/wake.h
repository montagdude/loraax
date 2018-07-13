// Header for Wake class

#ifndef WAKE_H
#define WAKE_H

#include <vector>
#include "vertex.h"
#include "tripanel.h"
#include "quadpanel.h"

class Panel;

/******************************************************************************/
//
// Wake class. Stores vortex rings and horseshoe vortices.
//
/******************************************************************************/
class Wake {

	private:

	int _idx;
	int _nspan, _nstream;						// Spanwise and streamwise
												//   vertices
	std::vector<Vertex> _verts;					// Vertices
	std::vector<double> _newx, _newy, _newz;	// New vertex positions after
	                                        	//   wake rollup
	std::vector<double> _newxinc, _newyinc, _newzinc;	
	std::vector<Vertex *> _topteverts, _botteverts;
												// Pointers to top and bottom TE
												//   vertices on wing
	std::vector<TriPanel> _tris;				// Tri doublet panels
	std::vector<QuadPanel> _quads;				// Quad doublet panels (trailing
												//   to near-infinity)

	public:

	// Constructor
	
	Wake ();
	
	// Initialize with upper and lower TE verts from the wing
	
	void initialize ( const std::vector<Vertex *> & topteverts,
	                  const std::vector<Vertex *> & botteverts,
	                  int & next_global_vertidx, int & next_global_elemidx,
	                  int wakeidx );

	// Get wake idx

	int idx () const;
	
	// Compute wake rollup and convect doublets downstream
	
	void convectVertices ( const double & dt,
	                       const std::vector<Panel *> & allsurf,
	                       const std::vector<Panel *> & allwake );
	void update ();

	// Access vertices and panels
	
	unsigned int nVerts () const;
	unsigned int nTris () const;
	unsigned int nQuads () const;
	Vertex * vert ( unsigned int vidx );
	TriPanel * triPanel ( unsigned int tidx );
	QuadPanel * quadPanel ( unsigned int qidx );

	// Induced velocity at a point due to planar wake aligned with freestream

	Eigen::Vector3d planarInducedVelocity ( const double & x, const double & y,
	                                        const double & z,
	                                        bool include_bound_leg=true,
	                                        int on_trailing_leg=-1) const;

	// Trefftz plane force calculation

	void farfieldForces ( const double & xtrefftz, const double & ztrefftz,
	                      const std::vector<Wake *> & allwake, double & lift,
	                      double & drag ) const;
};

#endif
