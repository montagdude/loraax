#include <vector>
#include <Eigen/Core>
#include <cmath>
#include "util.h"
#include "panel.h"
#include "vertex.h"

/******************************************************************************/
//
// Vertex class. Defines x, y, z coordinates in space and stores references to
// neighboring panels.
// 
/******************************************************************************/

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
Vertex::Vertex ()
{
  unsigned int i;

  _idx = -1;
  _x = 0.0;
  _y = 0.0;
  _z = 0.0;
  _npanels = 0;
  _panels.resize(0);
  _data.resize(8);
  for ( i = 0; i < 8; i++ ) { _data[i] = 0.; }
}

/******************************************************************************/
//
// Set/access index
//
/******************************************************************************/
void Vertex::setIdx ( int idx ) { _idx = idx; }
int Vertex::idx () const { return _idx; }

/******************************************************************************/
//
// Set or access coordinates
//
/******************************************************************************/
void Vertex::setCoordinates ( const double & x, const double & y,
                              const double & z )
{
  _x = x;
  _y = y;
  _z = z;
}

const double & Vertex::x () const { return _x; }
const double & Vertex::y () const { return _y; }
const double & Vertex::z () const { return _z; }

/******************************************************************************/
//
// Transformations
//
/******************************************************************************/
void Vertex::scale ( const double & factor )
{
  _x *= factor;
  _y *= factor;
  _z *= factor;
}

void Vertex::translate ( const double & dx, const double & dy,
                         const double & dz )
{
  _x += dx;
  _y += dy;
  _z += dz;
}

void Vertex::rotate ( const Eigen::Matrix3d & transform )
{
  Eigen::Vector3d point, transpoint;

  point(0) = _x;
  point(1) = _y;
  point(2) = _z;
  transpoint = transform*point;
  _x = transpoint(0);
  _y = transpoint(1);
  _z = transpoint(2);
}

/******************************************************************************/
//
// Adding or accessing panel references
//
/******************************************************************************/
int Vertex::addPanel ( Panel * panel )
{
  unsigned int i;

  // Only add if the panel is not already in the list

  for ( i = 0; i < _npanels; i++ )
  {
    if (_panels[i]->idx() == panel->idx())
    {
#ifdef DEBUG
      print_warning("Vertex::addPanel", "Panel " + 
                    int2string(panel->idx()) +
                    std::string(" already in list."));
#endif
      return 1;
    }
  }

  _npanels += 1;
  _panels.push_back(panel);

  return 0;
}

Panel * Vertex::panel ( unsigned int pidx )
{
  if (pidx >= _npanels)
  {
    conditional_stop(1, "Vertex::panel", "Index out of range.");
  }

  return _panels[pidx];
}

bool Vertex::isNeighbor ( const Panel * panel ) const
{
  unsigned int i;

  for ( i = 0; i < _npanels; i++ )
  {
    if (_panels[i]->idx() == panel->idx())
      return true;
  }

  return false;
}

unsigned int Vertex::nPanels () const { return _npanels; }

/******************************************************************************/
//
// Setting or accessing data. See legend in comments in vertex.h.
//
/******************************************************************************/
int Vertex::setData ( unsigned int idx, const double & var )
{
#ifdef DEBUG
  if (idx >= _data.size())
  {
    conditional_stop(1, "Vertex::setData", "Index out of range.");
    return 1;
  }
#endif

  _data[idx] = var;

  return 0;
}

const double & Vertex::data ( unsigned int idx ) const
{
#ifdef DEBUG
  if (idx >= _data.size())
    conditional_stop(1, "Vertex::data", "Index out of range.");
#endif

  return _data[idx];
}

/******************************************************************************/
//
// Computes data for surface vertices. Source strength, doublet strength, and
// velocity are interpolated from panels. Pressure and cp are computed from
// velocity using Bernoulli equation and compressibility corrections.
//
/******************************************************************************/
void Vertex::computeSurfaceData ()
{
  unsigned int i;
  double dx, dy, dz, dist, weightsum;
  Eigen::Vector3d cen;

  // Compute source strength, doublet strength, and velocity by interpolation
  // from panels

  _data[0] = 0.;
  _data[1] = 0.;
  _data[3] = 0.;
  _data[4] = 0.;
  _data[5] = 0.;
  weightsum = 0.;
  for ( i = 0; i < _npanels; i++ )
  {
    cen = _panels[i]->centroid();
    dx = _x - cen(0);
    dy = _y - cen(1);
    dz = _z - cen(2);
    dist = std::sqrt(std::pow(dx,2.) + std::pow(dy,2.) + std::pow(dz,2.));
    _data[0] += _panels[i]->sourceStrength()/dist;
    _data[1] += _panels[i]->doubletStrength()/dist;
    _data[3] += _panels[i]->velocity()(0)/dist;
    _data[4] += _panels[i]->velocity()(1)/dist;
    _data[5] += _panels[i]->velocity()(2)/dist;
    weightsum += 1./dist;
  }
  _data[0] /= weightsum;
  _data[1] /= weightsum;
  _data[3] /= weightsum;
  _data[4] /= weightsum;
  _data[5] /= weightsum;
}
