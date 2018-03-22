// Header for airfoil class

#ifndef AIRFOIL_H
#define AIRFOIL_H

#include <vector>
#include <string>
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
class Airfoil {

  private:

    // Geometric data

    unsigned int _nb, _n;		// Number of points (buffer, smoothed)
    std::vector<double> _xb, _zb;	// Input coordinates
    std::vector<double> _s, _xs, _zs;   // Spline fit of buffer coordinates
    double _sle;                        // Leading edge spline parameter
    std::vector<double> _x, _z;		// Smoothed coordinates
    double _y;				// y coordinate, used for sorting and
					//   interpolation

    bool _unit_transform;		// Transform flag (call unitTransform)

    // Aero data
    
    double _cl, _cd, _cm;
    std::vector<double> _cp, _cf;

  public:

    // Constructor

    Airfoil (); 

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
    bool splined () const;		// Whether spline fit has been done
    int splineInterp ( const double & sc, double & xc, double & zc ) const;
					// Coordinate at spline parameter

    // Scale buffer coordinates to unit chord and move to origin

    int unitTransform ();

    // Generate smoothed coordinates from buffer coordinates

    int smoothPaneling ( const xfoil_geom_options_type & geom_opts );

    // Airfoil data

    unsigned int nBuffer () const;
    unsigned int nSmoothed () const;
    void bufferCoordinates ( std::vector<double> & xb,
                             std::vector<double> & zb ) const;
    void smoothedCoordinates ( std::vector<double> & x,
                               std::vector<double> & z ) const;

    // Set/access y coordinate

    void setY ( const double & y );
    const double & y () const;
};

#endif
