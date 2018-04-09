#include <vector>
#include <Eigen/Core>
#include "util.h"
#include "vertex.h"
#include "singularities.h"
#include "vortex_ring.h"

/******************************************************************************/
//
// Vortex ring class
//
/******************************************************************************/

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
VortexRing::VortexRing ()
{
  _verts.resize(4);
  _type = "vortexring";
}

/******************************************************************************/
//
// Set vertices. Vertices should be given in a clockwise order.
//
/******************************************************************************/
int VortexRing::addVertex ( Vertex * vert )
{
  if (_currverts == 4)
  {
    conditional_stop(1, "VortexRing::addVertex",
                     "4 vertices have already been set.");
    return 1;
  }

  _verts[_currverts] = vert;
  _verts[_currverts]->addElement(this);
  _currverts += 1;

  return 0;
}

/******************************************************************************/
//
// Computes induced velocity at a point due to wake line with unit circulation
//
/******************************************************************************/
Eigen::Vector3d VortexRing::VCoeff ( const double & x, const double & y, 
                                     const double & z, const double & rcore,
                                     bool mirror_y ) const
{
  unsigned int i;
  Eigen::Vector3d vel, vel_mirror;

#ifdef DEBUG
  if (_currverts < 4)
    conditional_stop(1, "VortexRing::VCoeff",
                     "4 vertices have not yet been set.");
#endif

  // Initialize induced velocity

  vel(0) = 0.0;
  vel(1) = 0.0;
  vel(2) = 0.0;

  for ( i = 0; i < 3; i++ )
  {
    vel += vortex_velocity(x, y, z, _verts[i]->x(), _verts[i]->y(),
                           _verts[i]->z(), _verts[i+1]->x(), _verts[i+1]->y(),
                           _verts[i+1]->z(), rcore, false);
  }
  vel += vortex_velocity(x, y, z, _verts[3]->x(), _verts[3]->y(),
                         _verts[3]->z(), _verts[0]->x(), _verts[0]->y(),
                         _verts[0]->z(), rcore, false);

  // Compute mirror image velocity if requested

  if (mirror_y)
  {
    vel_mirror(0) = 0.0;
    vel_mirror(1) = 0.0;
    vel_mirror(2) = 0.0;

    for ( i = 0; i < 3; i++ )
    {
      vel_mirror += vortex_velocity(x, -y, z, _verts[i]->x(), _verts[i]->y(),
                                    _verts[i]->z(), _verts[i+1]->x(),
                                    _verts[i+1]->y(), _verts[i+1]->z(), rcore,
                                    false);
    }
    vel_mirror += vortex_velocity(x, -y, z, _verts[3]->x(), _verts[3]->y(),
                                  _verts[3]->z(), _verts[0]->x(),
                                  _verts[0]->y(), _verts[0]->z(), rcore, false);

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
Eigen::Vector3d VortexRing::inducedVelocity ( const double & x,
                                             const double & y, const double & z,
                                             const double & rcore,
                                             bool mirror_y ) const
{
  Eigen::Vector3d vel;

  vel = VCoeff(x, y, z, rcore, mirror_y)*_gamma; 
  return vel;
}
