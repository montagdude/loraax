// Header for aircraft class

#ifndef AIRCRAFT_H
#define AIRCRAFT_H

#include <vector>
#include <string>
#include <Eigen/Dense>
#include <fstream>
#include "wing.h"

class Vertex;
class Panel;

/******************************************************************************/
//
// Aircraft class. Contains some number of wings and related data and members.
//
/******************************************************************************/
class Aircraft {

  private:

    std::vector<Wing> _wings;		// Wings

    double _sref;			// Reference area
    double _lref;			// Pitching moment reference length
    Eigen::Vector3d _momcen;		// Pitching moment reference point

    std::vector<Vertex *> _verts;	// Pointers to vertices
    std::vector<Panel *> _panels;	// Pointers to panels 
    std::vector<Vertex *> _wakeverts;   // Pointers to wake vertices
    std::vector<Panel *> _wakepanels;	// Pointers to wake doublet panels

    Eigen::MatrixXd _sourceic, _doubletic;
					// Aero influence coefficients due to
					//   sources and doublets on surface
    Eigen::MatrixXd _aic;		// Aero influence coefficients matrix
    Eigen::VectorXd _mun;		// Normalized doublet strengths vector
    Eigen::VectorXd _rhs;		// Right hand side vector
    Eigen::PartialPivLU<Eigen::MatrixXd> _lu;
					// LU factorization of AIC matrix
    double _lift, _drag, _moment;	// Dimensional forces and moments
    double _cl, _cd, _cm;		// Force and moment coefficients

    // Set up pointers to vertices, panels, and wake elements

    void setGeometryPointers ();

    // Write VTK viz

    int writeSurfaceViz ( const std::string & fname ) const;
    void writeSurfaceData ( std::ofstream & f ) const;
    int writeWakeViz ( const std::string & fname ) const;
    void writeWakeData ( std::ofstream & f ) const;
    int writeWakeStripViz ( const std::string & prefix );

  public:

    // Constructor

    Aircraft ();

    // Read from XML

    int readXML ( const std::string & geom_file );

    // Set source and doublet strengths

    void setSourceStrengths ();
    void setDoubletStrengths ();
    void setWakeDoubletStrengths ( bool init );

    // Construct AIC matrix and RHS vector, factorize, and solve

    void constructSystem ( unsigned int iter );
    void factorize ();
    void solveSystem ();

    // Gives size of system of equations (= number of panels)

    unsigned int systemSize () const;

    // Computes surface velocities and pressures

    void computeSurfaceQuantities ();

    // Computes or access forces and moments

    void computePressureForceMoment ();
    const double & lift () const;
    const double & drag () const;
    const double & pitchingMoment () const;
    const double & liftCoefficient () const;
    const double & dragCoefficient () const;
    const double & pitchingMomentCoefficient () const;

    // Write forces and moments to file

    int writeForceMoment ( int iter ) const;

    // Computes sectional pressure forces

    void computeSectionPressureForces (); 

    // BL calculations with Xfoil

    void computeBL (); 

    // Write section forces to file

    void writeSectionForces ( int iter ) const;

    // Convects and updates wake panels

    void moveWake ();

    // Write VTK viz

    int writeViz ( const std::string & prefix, int iter ) const;
};

#endif
