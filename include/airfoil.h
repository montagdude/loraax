// Header for airfoil class

#ifndef AIRFOIL_H
#define AIRFOIL_H

#include <vector>
#include <string>
#include "sectional_object.h"
extern "C"
{
	#include <xfoil_interface.h>
}

/******************************************************************************/
//
// Airfoil class. Stores and manipulates airfoil coordinates and solves BL
// equations with Xfoil.
//
/******************************************************************************/
class Airfoil: public SectionalObject {

	private:

	// Xfoil data
	
	xfoil_data_group _xdg;
	
	// Geometric data
	
	int _nb, _n;						// Number of points (buffer, smoothed)
	std::vector<double> _s, _xs, _zs;	// Spline fit of buffer coordinates
	std::vector<double> _ssmoothed;		// Spline vector for smoothed airfoil
	double _sle;						// Leading edge spline parameter
	
	bool _unit_transform;				// Transform flag (call unitTransform)

	// Helper function to avoid reusing code between copy constructor and copy
	// assignment
	
	void copyData ( const Airfoil & foil );
	
	public:
	
	// Constructor, copy constructor, copy assignment, destructor
	
	Airfoil ();
	Airfoil ( const Airfoil & foil );
	Airfoil & operator = ( const Airfoil & foil );
	~Airfoil ();
	
	// Setting airfoil coordinates
	
	int readCoordinates ( const std::string & fname );
	int naca4Coordinates ( const std::string & des, const int & npointside );
	int naca5Coordinates ( const std::string & des, const int & npointside );
	int setCoordinates ( const std::vector<double> & x,
	                     const std::vector<double> & z );
	int interpCoordinates ( const Airfoil & foil1, const Airfoil & foil2,
	                        const double interpfrac );
	
	// Puts buffer coordinates in counterclockwise order
	
	void ccwOrderCoordinates ();
	
	// Spline fit coordinates
	
	int splineFit ();
	const double & sLen () const;	// Total arc length
	const double & sLE () const;	// LE spline coordinate
	bool splined () const;			// Whether spline fit has been done
	int splineInterp ( const double & sc, double & xc, double & zc ) const;
									// Coordinate at spline parameter
	
	// Scale buffer coordinates to unit chord and move to origin
	
	int unitTransform ();
	
	// Set xfoil options. Should do this before most other operations.
	
	void setXfoilOptions ( const xfoil_options_type & xfoil_opts,
	                       const xfoil_geom_options_type & geom_opts );
	
	// Generate smoothed coordinates from buffer coordinates
	
	int smoothPaneling ();
	
	// Smoothed s vector at given idx
	
	const double & sSmoothed ( int idx ) const;
	
	// Returns or modifies trailing edge gap. Note: modifyTEGap only only
	// modifies the buffer airfoil coordinates. You may want to call
	// splineInterp, smoothPaneling, etc. again afterwards.
	
	double teGap () const;
	int modifyTEGap ( const double & newgap, const double & blendloc );
	
	// Airfoil data
	
	int nBuffer () const;
	int nSmoothed () const;
	void bufferCoordinates ( std::vector<double> & xb,
	                         std::vector<double> & zb ) const;
	void smoothedCoordinates ( std::vector<double> & x,
	                           std::vector<double> & z ) const;

	// Running Xfoil and returning BL data
	
	void setReynoldsNumber ( const double & re );
	void setMachNumber ( const double & mach );
	int runXfoil ( const double & clspec );
	std::vector<double> blData ( const std::string & varname,
	                             int & stat ) const;

	// Wake data

	int nWake () const;
	void wakeCoordinates ( int nw, std::vector<double> & xw,
	                       std::vector<double> & zw ) const;
	std::vector<double> wakeSVector ( int nw ) const;
	std::vector<double> wakeDeltastar ( int nw ) const;
	std::vector<double> wakeUedge ( int nw ) const;
};

#endif
