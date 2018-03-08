// Header for airfoil class

#ifndef AIRFOIL_H
#define AIRFOIL_H

#include <vector>
#include <string>

/******************************************************************************/
//
// Airfoil class. Defines wing slice where BL equations are solved.
//
/******************************************************************************/
class Airfoil {

  private:

    unsigned int _nb, _n;		// Number of points (buffer, smoothed)
    double _chord;
    std::vector<double> _xb, _zb;	// Input coordinates
    std::vector<double> _x, _z;		// Smoothed coordinates

  public:

    // Constructor

    Airfoil (); 

    // Setting airfoil coordinates

    int readCoordinates ( const std::string & fname );
    int naca4Coordinates ( char des[4] ); 
    int naca5Coordinates ( char des[5] );
    int setCoordinates ( const std::vector<double> & x,
                         const std::vector<double> & z );
    int interpCoordinates ( const Airfoil & foil1, const Airfoil & foil2,
                            double interpfrac );

    // Puts buffer coordinates in counterclockwise order

    void ccOrderCoordinates ();
};

#endif
