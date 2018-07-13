// Header for Vertex class

#ifndef VERTEX_H
#define VERTEX_H

#include <vector>
#include <Eigen/Core>

class Panel;

/******************************************************************************/
//
// Vertex class. Defines x, y, z coordinates in space and stores references to
// neighboring panels.
//
/******************************************************************************/
class Vertex {

	private:

	int _idx;
	double _x, _y, _z;
	std::vector<Panel *> _panels;	// Neighboring panels
	unsigned int _npanels;
	bool _vizcoords;				// Separate coordinates for
	double _xviz, _yviz, _zviz;	 	//   visualization. Used for trailing
									//   wake vertices that are far
									//   downstream.
	bool _inccoords;				// Prandtl-Glauert transformed coordinates
	double _xinc, _yinc, _zinc;

	// Vertex data: source strength, doublet strength, Vx, Vy, Vz, pressure, cp,
	// cf, deltastar, ampl, uedge, cp from Xfoil
	
	std::vector<double> _data;
	
	double _waketime;			// Variable to track wake convection

	public:

	const static int dataSize = 12;     // Total number of data vars
	const static int firstBLData = 7;   // Index of first data var from BL calcs
	
	// Constructor
	
	Vertex ();
	
	// Set or access index
	
	void setIdx ( int idx );
	int idx () const;
	
	// Setting or accessing coordinates
	
	void setCoordinates ( const double & x, const double & y,
	                      const double & z );
	const double & x () const;
	const double & y () const;
	const double & z () const;
	
	void setVizCoordinates ( const double & x, const double & y,
	                         const double & z );
	const double & xViz () const;
	const double & yViz () const;
	const double & zViz () const;

	void setIncompressibleCoordinates ( const double & x, const double & y,
	                                    const double & z );
	const double & xInc () const;
	const double & yInc () const;
	const double & zInc () const;
	
	// Transformations
	
	void scale ( const double & factor );
	void translate ( const double & dx, const double & dy, const double & dz );
	void rotate ( const Eigen::Matrix3d & transform );
	
	// Adding or accessing connected panels
	
	int addPanel ( Panel * panel ); 
	Panel * panel ( unsigned int pidx );
	bool isNeighbor ( const Panel * panel ) const;
	unsigned int nPanels () const;
	
	// Setting or accessing data. See legend in comments above.
	
	int setData ( unsigned int idx, const double & var );
	const double & data ( unsigned int idx ) const;
	
	// Setting or accessing wake time variables
	
	void setWakeTime ( const double & waketime );
	void incrementWakeTime ( const double & dt );
	const double & wakeTime () const;
	
	// Averages neighboring panel quantities to vertex
	
	void averageFromPanels ();
	
	// Distance and vector to another vertex
	
	double distance ( const Vertex & vert ) const;
	Eigen::Vector3d vectorTo ( const Vertex & vertex ) const;
};

#endif
