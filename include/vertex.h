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
    bool _vizcoords;		 	// Separate coordinates for
    double _xviz, _yviz, _zviz;	 	//   visualization. Used for trailing
					//   wake vertices that are far
					//   downstream.

    // Vertex data: source strength, doublet strength, Vx, Vy, Vz, pressure, cp

    std::vector<double> _data;

  public:

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

    // Interpolates panel quantities to vertex

    void interpFromPanels ();

    // Distance and vector to another vertex

    double distance ( const Vertex & vert ) const;
    Eigen::Vector3d vectorTo ( const Vertex & vertex ) const;
};

#endif
