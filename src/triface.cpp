#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <Eigen/Core>
#include "util.h"
#include "settings.h"
#include "geometry.h"
#include "singularities.h"
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
  _x.resize(3);
  _y.resize(3);
  _z.resize(3);
  _xtrans.resize(3);
  _ytrans.resize(3);
}

/******************************************************************************/
//
// Sets endpoints and calculates face geometric quantities. Endpoints should
// be given in a counterclockwise (right-handed) direction.
//
/******************************************************************************/
void TriFace::setEndpoints ( const std::vector<double> & x,
                             const std::vector<double> & y,
                             const std::vector<double> & z )
{
  unsigned int i;

  if ( (x.size() != 3) || (y.size() != 3) || (z.size() != 3) )
  {
    conditional_stop(1, "TriFace::setEndpoints",
                     "x, y, and z must have 3 elements.");
  }

  for ( i = 0; i < 3; i++ )
  {
    _x[i] = x[i];
    _y[i] = y[i];
    _z[i] = z[i];
  }
  computeCharacteristicLength();
  computeArea();
  computeNormal();
  computeCentroid();
  computeTransform();
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

  side1(0) = _x[1] - _x[0];
  side1(1) = _y[1] - _y[0];
  side1(2) = _z[1] - _z[0];
  len1 = side1.norm();

  side2(0) = _x[2] - _x[1];
  side2(1) = _y[2] - _y[1];
  side2(2) = _z[2] - _z[1];
  len2 = side2.norm();

  side3(0) = _x[0] - _x[2];
  side3(1) = _y[0] - _y[2];
  side3(2) = _z[0] - _z[2];
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
  _area = tri_area(_x[0], _y[0], _z[0],
                   _x[1], _y[1], _z[1],
                   _x[2], _y[2], _z[2]);
}

/******************************************************************************/
//
// Computes normal
//
/******************************************************************************/
void TriFace::computeNormal ()
{
  _norm = tri_normal(_x[0], _y[0], _z[0],
                     _x[1], _y[1], _z[1],
                     _x[2], _y[2], _z[2]);
}

/******************************************************************************/
//
// Computes centroid by intersection of vertex - edge midpoint lines
//
/******************************************************************************/
void TriFace::computeCentroid () 
{
  _cen = tri_centroid(_x[0], _y[0], _z[0],
                      _x[1], _y[1], _z[1],
                      _x[2], _y[2], _z[2]);
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
    vec(0) = _x[i] - _cen[0];
    vec(1) = _y[i] - _cen[1];
    vec(2) = _z[i] - _cen[2];
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
