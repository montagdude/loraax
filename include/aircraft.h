// Header for aircraft class

#ifndef AIRCRAFT_H
#define AIRCRAFT_H

#include <vector>
#include <string>
#include <Eigen/Dense>
#include <fstream>
#include "wing.h"
#include "farfield.h"

class Vertex;
class Panel;
class Wake;

/******************************************************************************/
//
// Aircraft class. Contains some number of wings and related data and members.
//
/******************************************************************************/
class Aircraft {

    private:

    std::vector<Wing> _wings;           // Wings
    
    double _sref;                       // Reference area
    double _lref;                       // Pitching moment reference length
    Eigen::Vector3d _momcen;            // Pitching moment reference point

    double _xte, _zte, _maxspan;        // Furthest aft root TE points and
                                        //   largest wingspan, used for Trefftz
                                        //   plane calculation
    
    std::vector<Vertex *> _verts;       // Pointers to vertices
    std::vector<Panel *> _panels;       // Pointers to panels 
    std::vector<Vertex *> _wakeverts;   // Pointers to wake vertices
    std::vector<Panel *> _wakepanels;   // Pointers to wake doublet panels
    std::vector<Wake *> _allwake;       // Pointers to wakes

    Farfield _farfield;                 // Farfield (for post calculations only)
    
    Eigen::MatrixXd _sourceic, _doubletic;
                                        // Aero influence coefficients due to
                                        //   sources and doublets on surface
    Eigen::MatrixXd _aic;               // Aero influence coefficients matrix
    Eigen::VectorXd _mun;               // Normalized doublet strengths vector
    Eigen::VectorXd _rhs;               // Right hand side vector
    Eigen::PartialPivLU<Eigen::MatrixXd> _lu;
                                        // LU factorization of AIC matrix
    
    // Set up pointers to vertices, panels, and wake elements
    
    void setGeometryPointers ();
    
    // Write VTK viz
    
    int writeSurfaceViz ( const std::string & fname ) const;
    void writeSurfaceData ( std::ofstream & f ) const;
    void writeSurfaceScalar ( std::ofstream & f, const std::string & varname,
                              unsigned int varidx ) const;
    int writeWakeViz ( const std::string & fname ) const;
    void writeWakeData ( std::ofstream & f ) const;
    int writeWakeStripViz ( const std::string & prefix );

    public:

    // Constructor
    
    Aircraft ();
    
    // Read from XML
    
    int readXML ( const std::string & geom_file );
    
    // Set source and doublet strengths
    
    void setSourceStrengths ( bool init );
    void setDoubletStrengths ();
    
    // Construct AIC matrix and RHS vector, factorize, and solve
    
    void constructSystem ( bool init );
    void factorize ();
    void solveSystem ();
    
    // Gives size of system of equations (= number of panels)
    
    unsigned int systemSize () const;
    
    // Computes surface velocities and pressures
    
    void computeSurfaceQuantities ();
    
    // BL calculations with Xfoil
    
    void computeBL (); 

    // Sets up viscous wake for each wing

    void setupViscousWake ();
    
    // Compute or access forces and moments
    
    void computeForceMoment ();

    double lift () const;                       // Trefftz + skin friction
    double trefftzLift () const;                // Calculated in Trefftz plane
    double skinFrictionLift () const;           // Via skin friction integration
    double integratedLift () const;             // Via pressure integration

    double drag () const;                       // Induced + viscous
    double inducedDrag () const;                // Calculated in Trefftz plane
    double parasiticDrag () const;              // Via BL drag span integration
    double skinFrictionDrag () const;           // Skin fric. part of viscous
    double integratedDrag () const;             // Via pressure integration;
                                                //   inaccurate in viscous cases

    double pitchingMoment () const;             // Pressure + skin friction
    double pressurePitchingMoment () const;
    double skinFrictionPitchingMoment () const;

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
    
    // Write section force and moment coefficients to file
    
    void writeSectionForceMoment ( int iter ) const;
    
    // Convects and updates wake panels
    
    void moveWake ();
    
    // Write VTK viz
    
    int writeViz ( const std::string & prefix, int iter ) const;
};

#endif
