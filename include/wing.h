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
	double _splanform, _cbar, _span;// Planform area, MAC, and span
	
	std::vector<Section> _sections;	// Spanwise sections used for
									//   calculations, set by _nspan and
									//   root & tip spacing rations
	std::vector<double> _stations;	// Section positions in span coordinates
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

	double _liftp, _liftf;			// Pressure and skin friction, integrated
	double _lifttr;					// Trefftz plane
	double _dragp, _dragf;			// Pressure and skin friction, integrated
	double _dragv;					// Parasitic drag (from BL solution)
	double _dragtr;					// Trefftz plane
	double _momentp, _momentf;		// Pressure and skin friction, integrated
	
	std::vector<double> adjustSpacing ( 
	                           const std::vector<double> & nom_stations ) const;

	void computeAreaMAC ( const std::vector<Section> & sorted_user_sections );

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

	// Root trailing edge positions

	const double & xTE ();
	const double & zTE ();

	// Geometric quantities

	const double & span () const;
	const double & meanAerodynamicChord () const;
	const double & planformArea () const;
	
	// Creates panels and surface vertex pointers
	
	void createPanels ( int & next_global_vertidx, int & next_global_elemid );
	
	// Set up wake
	
	void setupWake ( const double & maxspan, int & next_global_vertidx,
	                 int & next_global_elemidx, int wakeidx );

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
	
	void computeForceMoment ( const Eigen::Vector3d & momcen,
	                          const double & xtrefftz, const double & ztrefftz,
	                          const std::vector<Wake *> & allwake );
	
	double lift () const;						// Trefftz + skin friction
	const double & trefftzLift () const;		// Calculated in Trefftz plane
	const double & skinFrictionLift () const;	// Via skin friction integration
	const double & integratedLift () const;		// Via pressure integration

	double drag () const;						// Induced + viscous
	const double & inducedDrag () const;		// Calculated in Trefftz plane
	const double & parasiticDrag () const;		// Via BL drag span integration
	const double & skinFrictionDrag () const;	// Skin fric. part of viscous
	const double & integratedDrag () const;		// Via pressure integration;
												//   inaccurate in viscous cases

	double pitchingMoment () const;				// Pressure + skin friction
	const double & pressurePitchingMoment () const;
	const double & skinFrictionPitchingMoment () const;

	double liftCoefficient () const;
	double trefftzLiftCoefficient () const;
	double skinFrictionLiftCoefficient () const;
	double integratedLiftCoefficient () const;

	double dragCoefficient () const;
	double inducedDragCoefficient () const;
	double parasiticDragCoefficient () const;
	double skinFrictionDragCoefficient () const;
	double integratedDragCoefficient () const;

	double pitchingMomentCoefficient () const;
	double pressurePitchingMomentCoefficient () const;
	double skinFrictionPitchingMomentCoefficient () const;
	
	// Write forces and moments to file
	
	int writeForceMoment ( int iter ) const;
	
	// Write sectional force and moment coefficients to file
	
	int writeSectionForceMoment ( int iter ) const;
};

#endif
