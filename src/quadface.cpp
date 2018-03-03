#include <vector>
#include <string>
#include <cmath>
#include <Eigen/Core>
#include "util.h"
#include "settings.h"
#include "geometry.h"
#include "singularities.h"
#include "node.h"
#include "edge.h"
#include "face.h"
#include "quadface.h"

/******************************************************************************/
//
// Quad face class.  Derived class from Face, has 4 nodes/edges.
// 
/******************************************************************************/

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
QuadFace::QuadFace ()
{
  _noderef.resize(4);
  _edgeref.resize(4);
  _xtrans.resize(4);
  _ytrans.resize(4);
}

/******************************************************************************/
//
// Computes face geometric quantities once four nodes are set
//
/******************************************************************************/
void QuadFace::setNode ( unsigned int nidx, Node * nodein )
{
  Face::setNode(nidx, nodein);
  if (_currnodes == 4)
  {
    computeCharacteristicLength();
    computeArea();
    computeNormal();
    computeCentroid();
    computeTransform();
  }
}

/******************************************************************************/
//
// Computes characteristic length as max diagaonal length (used for farfield
// source / doublet approximations)
//
/******************************************************************************/
void QuadFace::computeCharacteristicLength ()
{
  Eigen::Vector3d diag1, diag2;

#ifndef NDEBUG

  if (_currnodes < 4)
  {
    conditional_stop(1, "QuadFace::computeCharacteristicLength",
                     "Not enough nodes set for quadrilateral face.");
  }

#endif

  // Diagonals

  diag1(0) = node(2).x() - node(0).x();
  diag1(1) = node(2).y() - node(0).y();
  diag1(2) = node(2).z() - node(0).z();

  diag2(0) = node(3).x() - node(1).x();
  diag2(1) = node(3).y() - node(1).y();
  diag2(2) = node(3).z() - node(1).z();

  _length = std::max(diag1.norm(), diag2.norm());
}

/******************************************************************************/
//
// Computes area
//
/******************************************************************************/
void QuadFace::computeArea ()
{
#ifndef NDEBUG

  if (_currnodes < 4)
  {
    conditional_stop(1, "QuadFace::computeArea",
                     "Not enough nodes set for quadrilateral face.");
  }

#endif

  _area = quad_area(node(0).x(), node(0).y(), node(0).z(),
                    node(1).x(), node(1).y(), node(1).z(),
                    node(2).x(), node(2).y(), node(2).z(),
                    node(3).x(), node(3).y(), node(3).z());
}

/******************************************************************************/
//
// Computes normal
//
/******************************************************************************/
void QuadFace::computeNormal ()
{
#ifndef NDEBUG

  if (_currnodes < 4)
  {
    conditional_stop(1, "QuadFace::computeNormal",
                     "Not enough nodes set for quadrilateral face.");
  }

#endif

  _norm = quad_normal(node(0).x(), node(0).y(), node(0).z(),
                      node(1).x(), node(1).y(), node(1).z(),
                      node(2).x(), node(2).y(), node(2).z(),
                      node(3).x(), node(3).y(), node(3).z());
}

/******************************************************************************/
//
// Computes centroid by dividing into two triangles
//
/******************************************************************************/
void QuadFace::computeCentroid () 
{
#ifndef NDEBUG

  if (_currnodes < 4)
  {
    conditional_stop(1, "QuadFace::computeCentroid",
                     "Not enough nodes set for quadrilateral face.");
  }

#endif

  _cen = quad_centroid(node(0).x(), node(0).y(), node(0).z(),
                       node(1).x(), node(1).y(), node(1).z(),
                       node(2).x(), node(2).y(), node(2).z(),
                       node(3).x(), node(3).y(), node(3).z());
}

/******************************************************************************/
//
// Computes transform from inertial frame to panel frame and inverse transform.
// Also computes and stores panel endpoints in panel frame.
//
/******************************************************************************/
void QuadFace::computeTransform ()
{
  unsigned int i;
  Eigen::Vector3d vec, transvec;

  // Transformation from inertial frame to panel frame

  _trans = transform_from_normal(_norm[0], _norm[1], _norm[2]);

  // Inverse transform

  _invtrans = _trans.transpose();

  // Convert panel endpoints to panel frame (only store x and y, since z
  // should be 0 or nearly 0)

  for ( i = 0; i < 4; i++ )
  {
    vec(0) = node(i).x() - _cen[0];
    vec(1) = node(i).y() - _cen[1];
    vec(2) = node(i).z() - _cen[2];
    transvec = _trans*vec;
    _xtrans[i] = transvec(0);
    _ytrans[i] = transvec(1);
  }
}

/******************************************************************************/
//
// Computes source influence coefficient
//
/******************************************************************************/
double QuadFace::sourcePhiCoeff ( const double & x, const double & y, 
                                  const double & z, const bool onpanel,
                                  const std::string & side ) const
{
  Eigen::Vector3d vec, transvec;

#ifndef NDEBUG

  if (_currnodes < 4)
  {
    conditional_stop(1, "QuadFace::sourcePhiCoeff",
                     "Not enough nodes set for quadrilateral face.");
  }

#endif

  // Vector from centroid to point

  vec(0) = x - _cen[0];
  vec(1) = y - _cen[1];
  vec(2) = z - _cen[2];

  // Determine whether to use farfield approximation

  if (vec.norm() > farfield_distance_factor*_length)
  {
    // Transformation not needed for point source

    return _area*point_source_potential(vec(0), vec(1), vec(2));
  }

  else
  {
    // Transform point to panel frame
  
    transvec = _trans*vec;
  
    // Source potential -- panel endpoints given in clockwise order
  
    return quad_source_potential(transvec(0), transvec(1), transvec(2),
                                 _xtrans[0], _ytrans[0], _xtrans[3], _ytrans[3], 
                                 _xtrans[2], _ytrans[2], _xtrans[1], _ytrans[1],
                                 onpanel, side);
  }
}

/******************************************************************************/
//
// Computes source induced velocity component
//
/******************************************************************************/
Eigen::Vector3d QuadFace::sourceVCoeff ( const double & x, const double & y, 
                                         const double & z, const bool onpanel,
                                         const std::string & side ) const
{
  Eigen::Vector3d vec, transvec, velpf, velif;

#ifndef NDEBUG

  if (_currnodes < 4)
  {
    conditional_stop(1, "QuadFace::sourceVCoeff",
                     "Not enough nodes set for quadrilateral face.");
  }

#endif

  // Vector from centroid to point

  vec(0) = x - _cen[0];
  vec(1) = y - _cen[1];
  vec(2) = z - _cen[2];

  // Determine whether to use farfield approximation 

  if (vec.norm() > farfield_distance_factor*_length)
  {
    // Transformation not needed for point source

    velif = _area*point_source_velocity(vec(0), vec(1), vec(2));
  }

  else
  {
    // Transform point to panel frame
  
    transvec = _trans*vec;
  
    // Source velocity -- panel endpoints given in clockwise order
  
    velpf = quad_source_velocity(transvec(0), transvec(1), transvec(2),
                                 _xtrans[0], _ytrans[0], _xtrans[3], _ytrans[3], 
                                 _xtrans[2], _ytrans[2], _xtrans[1], _ytrans[1],
                                 onpanel, side);
  
    // Convert to inertial frame
  
    velif = _invtrans*velpf;
  }

  return velif;
}

/******************************************************************************/
//
// Computes doublet influence coefficient
//
/******************************************************************************/
double QuadFace::doubletPhiCoeff ( const double & x, const double & y, 
                                   const double & z, const bool onpanel,
                                   const std::string & side ) const
{
  Eigen::Vector3d vec, transvec;

#ifndef NDEBUG

  if (_currnodes < 4)
  {
    conditional_stop(1, "QuadFace::doubletPhiCoeff",
                     "Not enough nodes set for quadrilateral face.");
  }

#endif

  // Transform point to panel frame

  vec(0) = x - _cen[0];
  vec(1) = y - _cen[1];
  vec(2) = z - _cen[2];
  transvec = _trans*vec;

  // Determine whether to use farfield approximation

  if (vec.norm() > farfield_distance_factor*_length)
  {
    return _area*point_doublet_potential(transvec(0), transvec(1), transvec(2));
  }

  else
  {
    // Doublet potential -- points given in clockwise order
  
    return quad_doublet_potential(transvec(0), transvec(1), transvec(2),
                                 _xtrans[0], _ytrans[0], _xtrans[3], _ytrans[3],
                                 _xtrans[2], _ytrans[2], _xtrans[1], _ytrans[1],
                                 onpanel, side);
  }
}

/******************************************************************************/
//
// Computes doublet induced velocity component
//
/******************************************************************************/
Eigen::Vector3d QuadFace::doubletVCoeff ( const double & x, const double & y, 
                                          const double & z, const bool onpanel,
                                          const std::string & side ) const
{
  Eigen::Vector3d vec, transvec, velpf, velif;

#ifndef NDEBUG

  if (_currnodes < 4)
  {
    conditional_stop(1, "QuadFace::doubletVCoeff",
                     "Not enough nodes set for quadrilateral face.");
  }

#endif

  // Transform point to panel frame

  vec(0) = x - _cen[0];
  vec(1) = y - _cen[1];
  vec(2) = z - _cen[2];
  transvec = _trans*vec;

  // Determine whether to use farfield approximation

  if (vec.norm() > farfield_distance_factor*_length)
  {
    velpf = _area*point_doublet_velocity(transvec(0), transvec(1), transvec(2));
  }

  else
  {
    // Doublet velocity -- points given in clockwise order
  
    velpf = quad_doublet_velocity(transvec(0), transvec(1), transvec(2),
                                 _xtrans[0], _ytrans[0], _xtrans[3], _ytrans[3],
                                 _xtrans[2], _ytrans[2], _xtrans[1], _ytrans[1],
                                 onpanel, side);
  }

  // Convert to inertial frame

  velif = _invtrans*velpf;

  return velif;
}

/******************************************************************************/
//
// Computes induced velocity potential at a point
//
/******************************************************************************/
double QuadFace::inducedPotential ( const double & x, const double & y, 
                                    const double & z, const bool onpanel,
                                    const std::string & side ) const
{
  double potential;

  potential = ( sourcePhiCoeff(x, y, z, onpanel, side)*_sigma
            +   doubletPhiCoeff(x, y, z, onpanel, side)*_mu );
  
  return potential;
}

/******************************************************************************/
//
// Computes induced velocity component at a point
//
/******************************************************************************/
Eigen::Vector3d QuadFace::inducedVelocity ( 
                                           const double & x, const double & y, 
                                           const double & z, const bool onpanel,
                                           const std::string & side ) const
{
  Eigen::Vector3d vel;

  vel = ( sourceVCoeff(x, y, z, onpanel, side)*_sigma
      +   doubletVCoeff(x, y, z, onpanel, side)*_mu );

  return vel;
}

