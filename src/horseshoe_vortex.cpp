#include <vector>
#include <Eigen/Core>
#include "util.h"
#include "vertex.h"
#include "singularities.h"
#include "horseshoe_vortex.h"

/******************************************************************************/
//
// Horseshoe vortex class
//
/******************************************************************************/

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
HorseshoeVortex::HorseshoeVortex ()
{
  _verts.resize(4);
  _type = "horseshoevortex";
}

/******************************************************************************/
//
// Set vertices. Vertices should be given in a clockwise order.
//
/******************************************************************************/
int HorseshoeVortex::addVertex ( Vertex * vert, bool ref_element_to_vert )
{
  if (_currverts == 4)
  {
    conditional_stop(1, "HorseshoeVortex::addVertex",
                     "4 vertices have already been set.");
    return 1;
  }

  _verts[_currverts] = vert;
  _currverts += 1;

  return 0;
}

/******************************************************************************/
//
// Computes induced velocity at a point due to wake line with unit circulation
//
/******************************************************************************/
Eigen::Vector3d HorseshoeVortex::VCoeff ( const double & x, const double & y, 
                                         const double & z, const double & rcore,
                                         bool mirror_y ) const
{
  Eigen::Vector3d vel, vel_mirror;

#ifdef DEBUG
  if (_currverts < 4)
    conditional_stop(1, "HorseshoeVortex::VCoeff",
                     "4 vertices have not yet been set.");
#endif

  // Initialize induced velocity

  vel(0) = 0.0;
  vel(1) = 0.0;
  vel(2) = 0.0;

  // First leg given in reverse order and negated, since point 1 is at infinity

  vel -= vortex_velocity(x, y, z, _verts[1]->x(), _verts[1]->y(),
                         _verts[1]->z(), _verts[0]->x(), _verts[0]->y(),
                         _verts[0]->z(), rcore, true);
  vel += vortex_velocity(x, y, z, _verts[1]->x(), _verts[1]->y(),
                         _verts[1]->z(), _verts[2]->x(), _verts[2]->y(),
                         _verts[2]->z(), rcore, false);
  vel += vortex_velocity(x, y, z, _verts[2]->x(), _verts[2]->y(),
                         _verts[2]->z(), _verts[3]->x(), _verts[3]->y(),
                         _verts[3]->z(), rcore, true);

  if (mirror_y)
  {
    vel_mirror(0) = 0.0;
    vel_mirror(1) = 0.0;
    vel_mirror(2) = 0.0;

    vel_mirror -= vortex_velocity(x, -y, z, _verts[1]->x(), _verts[1]->y(),
                                  _verts[1]->z(), _verts[0]->x(),
                                  _verts[0]->y(), _verts[0]->z(), rcore, true);
    vel_mirror += vortex_velocity(x, -y, z, _verts[1]->x(), _verts[1]->y(),
                                  _verts[1]->z(), _verts[2]->x(),
                                  _verts[2]->y(), _verts[2]->z(), rcore, false);
    vel_mirror += vortex_velocity(x, -y, z, _verts[2]->x(), _verts[2]->y(),
                                  _verts[2]->z(), _verts[3]->x(),
                                  _verts[3]->y(), _verts[3]->z(), rcore, true);

    vel(0) += vel_mirror(0);
    vel(1) -= vel_mirror(1);
    vel(2) += vel_mirror(2);
  }

  return vel;
}

/******************************************************************************/
//
// Computes induced velocity at a point
//
/******************************************************************************/
Eigen::Vector3d HorseshoeVortex::inducedVelocity ( const double & x,
                                             const double & y, const double & z,
                                             const double & rcore,
                                             bool mirror_y ) const
{
  Eigen::Vector3d vel;

  vel = VCoeff(x, y, z, rcore, mirror_y)*_gamma; 
  return vel;
}
