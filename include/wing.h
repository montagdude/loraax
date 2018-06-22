// Header for wing class

#ifndef WING_H
#define WING_H

#include <vector>
#include <string>
#include "section.h"
#include "airfoil.h"
#include "vertex.h"
#include "panel.h"
#include "quadpanel.h"
#include "tripanel.h"
#include "wake.h"
#include "wake_strip.h"
#include "viscous_wake.h"

/******************************************************************************/
//
// Wing class. Contains sections, airfoils, panels, wake, etc.
//
/******************************************************************************/
class Wing {

	private:

	std::string _name;
	unsigned int _nchord, _nspan;	// Number of points along chord and
									//   span. nchord is on each side (top
									//   and bottom); total is 2*nchord-1.
	unsigned int _ntipcap;			// Number of points along the arc
									//   forming the tip cap. Now fixed at
									//   3, and tips are not rounded.
	double _lesprat, _tesprat;		// LE and TE spacing ratios relative to
									//   uniform
	double _rootsprat, _tipsprat;	// Root and tip spacing ratios relative
									//   to uniform
	
	std::vector<Section> _sections;	// Spanwise sections used for
									//   calculations, set by _nspan and
									//   root & tip spacing rations
	std::vector<Airfoil> _foils;  	// User-specified airfoils
	std::vector<Vertex *> _verts;	// Pointers to vertices on wing
	std::vector<std::vector<Vertex> > _tipverts;
									// Vertices on interior of wing cap
	std::vector<QuadPanel> _quads;	// Quad panels
	std::vector<TriPanel> _tris;	// Tri panels at tip
	std::vector<std::vector<Panel *> > _panels;
									// Pointers to panels arranged in i-j
									//   coordinate directions
	Wake _wake;						// Wake
	std::vector<WakeStrip> _wakestrips;	
									// Wake strips behind TE panels
	ViscousWake _vwake;				// Viscous wake

	double _liftp, _liftv;			// Dimensional forces and moments
	double _dragp, _dragv;
	double _momentp, _momentv;
	double _clp, _clv, _cdp, _cdv;	// Force coefficients
	double _cmp, _cmv;				// Moment coefficients
	
	std::vector<double> adjustSpacing ( 
	                           const std::vector<double> & nom_stations ) const;
	
	void computeAreaMAC ( const std::vector<Section> &
	                      sorted_user_sections ) const;

	public:

	// Constructor
	
	Wing ();
	
	// Set/get name
	
	void setName ( const std::string & name );
	const std::string & name () const;
	
	// Set discretization and spacing info
	
	void setDiscretization ( unsigned int nchord, unsigned int nspan,
	                         const double & lesprat, const double & tesprat,
	                         const double & rootsprat,
	                         const double & tipsprat );
	
	// Set airfoils. Section airfoils are interpolated from these airfoils.
	// Airfoils do not need to be given in sorted order.
	
	void setAirfoils ( std::vector<Airfoil> & foils );
	
	// Set sections based on user section inputs and spacing. User sections
	// define the planform shape, and the _sections variables defines the final
	// discretization. Sections do not need to be given in sorted order.
	// Also computes planform area and MAC and writes to stdout.
	
	int setupSections ( std::vector<Section> & user_sections );
	
	// Creates panels and surface vertex pointers
	
	void createPanels ( int & next_global_vertidx, int & next_global_elemid );
	
	// Set up wake
	
	void setupWake ( int & next_global_vertidx, int & next_global_elemidx );

	// Computes velocities and pressures on surface panels; interpolate to
	// vertices
	
	void computeSurfaceQuantities ();
	
	// Access to verts and panels
	
	unsigned int nVerts () const;
	unsigned int nQuads () const;
	unsigned int nTris () const;
	Vertex * vert ( unsigned int vidx );
	QuadPanel * quadPanel ( unsigned int qidx );
	TriPanel * triPanel ( unsigned int tidx );
	
	// Access to wake and wake strips
	
	Wake & wake ();
	unsigned int nWStrips () const;
	WakeStrip * wStrip ( unsigned int wsidx );
	ViscousWake & viscousWake ();

	// Compute viscous forces (and skin friction, etc.) using Xfoil at sections
	
	void computeBL ();
	
	// Set up viscous wake (note: only possible after computing BL first time)
	
	void setupViscousWake ( int & next_global_vertidx,
	                        int & next_global_elemidx );

	// Compute forces and moments, including sectional
	
	void computeForceMoment ( const double & sref, const double & lref,
	                          const Eigen::Vector3d & momcen );
	
	double lift () const;
	const double & pressureLift () const;
	const double & viscousLift () const;
	
	double drag () const;
	const double & pressureDrag () const;
	const double & viscousDrag () const;
	
	double pitchingMoment () const;
	const double & pressurePitchingMoment () const;
	const double & viscousPitchingMoment () const;
	
	double liftCoefficient () const;
	const double & pressureLiftCoefficient () const;
	const double & viscousLiftCoefficient () const;
	
	double dragCoefficient () const;
	const double & pressureDragCoefficient () const;
	const double & viscousDragCoefficient () const;
	
	double pitchingMomentCoefficient () const;
	const double & pressurePitchingMomentCoefficient () const;
	const double & viscousPitchingMomentCoefficient () const;
	
	// Write forces and moments to file
	
	int writeForceMoment ( int iter ) const;
	
	// Write sectional force and moment coefficients to file
	
	int writeSectionForceMoment ( int iter ) const;
};

#endif
