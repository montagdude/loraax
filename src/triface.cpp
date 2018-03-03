#include <iostream>
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
#include "triface.h"

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
TriFace::TriFace ()
{
  _noderef.resize(3);
  _edgeref.resize(3);
  _xtrans.resize(3);
  _ytrans.resize(3);
}

/******************************************************************************/
//
// Computes face geometric quantities once three nodes are set
//
/******************************************************************************/
void TriFace::setNode ( unsigned int nidx, Node * nodein )
{
  Face::setNode(nidx, nodein);
  if (_currnodes == 3)
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
// Computes characteristic length as max side length (used for farfield 
// source / doublet approximations
//
/******************************************************************************/
void TriFace::computeCharacteristicLength ()
{
  double len1, len2, len3;
  Eigen::Vector3d side1, side2, side3;

#ifndef NDEBUG

  if (_currnodes < 3)
  {
    conditional_stop(1, "TriFace::computeCharacteristicLength",
                     "Not enough nodes set for triangular face.");
  }

#endif

  // Side lengths

  side1(0) = node(1).x() - node(0).x();
  side1(1) = node(1).y() - node(0).y();
  side1(2) = node(1).z() - node(0).z();
  len1 = side1.norm();

  side2(0) = node(2).x() - node(1).x();
  side2(1) = node(2).y() - node(1).y();
  side2(2) = node(2).z() - node(1).z();
  len2 = side2.norm();

  side3(0) = node(0).x() - node(2).x();
  side3(1) = node(0).y() - node(2).y();
  side3(2) = node(0).z() - node(2).z();
  len3 = side3.norm();

  _length = std::max(len1, std::max(len2, len3));
}

/******************************************************************************/
//
// Computes area
//
/******************************************************************************/
void TriFace::computeArea ()
{
#ifndef NDEBUG

  if (_currnodes < 3)
  {
    conditional_stop(1, "TriFace::computeArea",
                     "Not enough nodes set for triangular face.");
  }

#endif

  _area = tri_area(node(0).x(), node(0).y(), node(0).z(),
                   node(1).x(), node(1).y(), node(1).z(),
                   node(2).x(), node(2).y(), node(2).z());
}

/******************************************************************************/
//
// Computes normal
//
/******************************************************************************/
void TriFace::computeNormal ()
{
#ifndef NDEBUG

  if (_currnodes < 3)
  {
    conditional_stop(1, "TriFace::computeNormal",
                     "Not enough nodes set for triangular face.");
  }

#endif

  _norm = tri_normal(node(0).x(), node(0).y(), node(0).z(),
                     node(1).x(), node(1).y(), node(1).z(),
                     node(2).x(), node(2).y(), node(2).z());
}

/******************************************************************************/
//
// Computes centroid by intersection of vertex - edge midpoint lines
//
/******************************************************************************/
void TriFace::computeCentroid () 
{
#ifndef NDEBUG

  if (_currnodes < 3)
  {
    conditional_stop(1, "TriFace::computeCentroid",
                     "Not enough nodes set for triangular face.");
  }

#endif

  _cen = tri_centroid(node(0).x(), node(0).y(), node(0).z(),
                      node(1).x(), node(1).y(), node(1).z(),
                      node(2).x(), node(2).y(), node(2).z());
}

/******************************************************************************/
//
// Computes transform from inertial frame to panel frame and inverse transform.
// Also computes and stores panel endpoints in panel frame.
//
/******************************************************************************/
void TriFace::computeTransform ()
{
  unsigned int i;
  Eigen::Vector3d vec, transvec;

  // Transformation from inertial frame to panel frame

  _trans = transform_from_normal(_norm[0], _norm[1], _norm[2]);

  // Inverse transform

  _invtrans = _trans.transpose();

  // Convert panel endpoints to panel frame (only store x and y, since z
  // should be 0 or nearly 0)

  for ( i = 0; i < 3; i++ )
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
double TriFace::sourcePhiCoeff ( const double & x, const double & y, 
                                 const double & z, const bool onpanel,
                                 const std::string & side ) const
{
  Eigen::Vector3d vec, transvec;

#ifndef NDEBUG

  if (_currnodes < 3)
  {
    conditional_stop(1, "TriFace::sourcePhiCoeff",
                     "Not enough nodes set for triangular face.");
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
  
    return tri_source_potential(transvec(0), transvec(1), transvec(2),
                                _xtrans[0], _ytrans[0], _xtrans[2], _ytrans[2], 
                                _xtrans[1], _ytrans[1], onpanel, side);
  }
}

/******************************************************************************/
//
// Computes source induced velocity component
//
/******************************************************************************/
Eigen::Vector3d TriFace::sourceVCoeff ( const double & x, const double & y, 
                                        const double & z, const bool onpanel,
                                        const std::string & side ) const
{
  Eigen::Vector3d vec, transvec, velpf, velif;

#ifndef NDEBUG

  if (_currnodes < 3)
  {
    conditional_stop(1, "TriFace::sourceVCoeff",
                     "Not enough nodes set for triangular face.");
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
  
    velpf = tri_source_velocity(transvec(0), transvec(1), transvec(2),
                                _xtrans[0], _ytrans[0], _xtrans[2], _ytrans[2], 
                                _xtrans[1], _ytrans[1], onpanel, side);
  
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
double TriFace::doubletPhiCoeff ( const double & x, const double & y, 
                                  const double & z, const bool onpanel,
                                  const std::string & side ) const
{
  Eigen::Vector3d vec, transvec;

#ifndef NDEBUG

  if (_currnodes < 3)
  {
    conditional_stop(1, "TriFace::doubletPhiCoeff",
                     "Not enough nodes set for triangular face.");
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
  
    return tri_doublet_potential(transvec(0), transvec(1), transvec(2),
                                 _xtrans[0], _ytrans[0], _xtrans[2], _ytrans[2],
                                 _xtrans[1], _ytrans[1], onpanel, side);
  }
}

/******************************************************************************/
//
// Computes doublet induced velocity component
//
/******************************************************************************/
Eigen::Vector3d TriFace::doubletVCoeff ( const double & x, const double & y, 
                                         const double & z, const bool onpanel,
                                         const std::string & side ) const
{
  Eigen::Vector3d vec, transvec, velpf, velif;

#ifndef NDEBUG

  if (_currnodes < 3)
  {
    conditional_stop(1, "TriFace::doubletVCoeff",
                     "Not enough nodes set for triangular face.");
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
  
    velpf = tri_doublet_velocity(transvec(0), transvec(1), transvec(2),
                                 _xtrans[0], _ytrans[0], _xtrans[2], _ytrans[2],
                                 _xtrans[1], _ytrans[1], onpanel, side);
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
double TriFace::inducedPotential ( const double & x, const double & y, 
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
Eigen::Vector3d TriFace::inducedVelocity ( const double & x, const double & y, 
                                           const double & z, const bool onpanel,
                                           const std::string & side ) const
{
  Eigen::Vector3d vel;

  vel = ( sourceVCoeff(x, y, z, onpanel, side)*_sigma
      +   doubletVCoeff(x, y, z, onpanel, side)*_mu );

  return vel;
}

