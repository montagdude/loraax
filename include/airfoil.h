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
// Airfoil class. Defines wing slice where BL equations are solved.
//
/******************************************************************************/
class Airfoil {

  private:

    // Geometric data

    unsigned int _nb, _n;		// Number of points (buffer, smoothed)
    std::vector<double> _xb, _zb;	// Input coordinates
    std::vector<double> _x, _z;		// Smoothed coordinates

    // Geometric scaling/transformation data

    bool _unit_scale;                   // Scaling flag
    double _chord, _xle, _zle, _pitch;

    // Aero data
    
    double _cl, _cd, _cm;
    std::vector<double> _cp, _cf;

  public:

    // Constructor

    Airfoil (); 

    // Setting airfoil coordinates

    int readCoordinates ( const std::string & fname );
    int naca4Coordinates ( const std::string & des, int & npointside );
    int naca5Coordinates ( const std::string & des, int & npointside );
    int setCoordinates ( const std::vector<double> & x,
                         const std::vector<double> & z );
    int interpCoordinates ( const Airfoil & foil1, const Airfoil & foil2,
                            const double interpfrac );

    // Puts buffer coordinates in counterclockwise order

    void ccwOrderCoordinates ();

    // Scale buffer coordinates to unit chord and move to origin

    int unitTransform ();

    // Generate smoothed coordinates from buffer coordinates

    int smoothPaneling ( const xfoil_geom_options_type & geom_opts );

    // Set / retrieve scaling information

    int setScaling ( const double & chord, const double & xle,
                     const double & zle, const double & pitch );
    const double & chord () const;
    const double & xLE () const;
    const double & zLE () const;
    const double & pitch () const;

    // Airfoil data

    unsigned int nBuffer () const;
    unsigned int nSmoothed () const;
    void bufferCoordinates ( std::vector<double> & xb,
                             std::vector<double> & zb ) const;
    void smoothedCoordinates ( std::vector<double> & x,
                               std::vector<double> & z ) const;
};

#endif
