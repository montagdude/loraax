#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <Eigen/Core>
#include "settings.h"
#include "util.h"
#include "geometry.h"
#include "singularities.h"
#include "vertex.h"
#include "panel.h"
#include "tripanel.h"

/******************************************************************************/
//
// Quad panel class.  Derived class from Panel, has 4 nodes/edges.
// 
/******************************************************************************/

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
TriPanel::TriPanel ()
{
  _verts.resize(3);
  _xtrans.resize(3);
  _ytrans.resize(3);
  _type = "tripanel";
}

/******************************************************************************/
//
// Add vertex and compute panel geometric quantities when three are set.
// Vertices should be added in a counterclockwise order.
//
/******************************************************************************/
int TriPanel::addVertex ( Vertex * vert, bool ref_element_to_vert )
{
  if (_currverts == 3)
  {
    conditional_stop(1, "TriPanel::addVertex",
                     "3 vertices have already been set.");
    return 1;
  }

  _verts[_currverts] = vert;
  if (ref_element_to_vert)
    _verts[_currverts]->addPanel(this);
  _currverts += 1;

  if (_currverts == 3)
    recomputeGeometry();

  return 0;
}

/******************************************************************************/
//
// Recomputes all geometric quantities
//
/******************************************************************************/
int TriPanel::recomputeGeometry ()
{
#ifdef DEBUG
  if (_currverts != 3)
  {
    conditional_stop(1, "TriPanel::recomputeGeometry",
                     "3 vertices are required.");
    return 1;
  }
#endif

  computeCharacteristicLength();
  computeArea();
  computeNormal();
  computeCentroid();
  computeTransform();

  return 0;
}

/******************************************************************************/
//
// Computes characteristic length as max side length (used for farfield 
// source / doublet approximations
//
/******************************************************************************/
void TriPanel::computeCharacteristicLength ()
{
  double len1, len2, len3;
  Eigen::Vector3d side1, side2, side3;

  // Side lengths

  side1(0) = _verts[1]->x() - _verts[0]->x();
  side1(1) = _verts[1]->y() - _verts[0]->y();
  side1(2) = _verts[1]->z() - _verts[0]->z();
  
  len1 = side1.norm();

  side2(0) = _verts[2]->x() - _verts[1]->x();
  side2(1) = _verts[2]->y() - _verts[1]->y();
  side2(2) = _verts[2]->z() - _verts[1]->z();
  len2 = side2.norm();

  side3(0) = _verts[0]->x() - _verts[2]->x();
  side3(1) = _verts[0]->y() - _verts[2]->y();
  side3(2) = _verts[0]->z() - _verts[2]->z();
  len3 = side3.norm();

  _length = std::max(len1, std::max(len2, len3));
}

/******************************************************************************/
//
// Computes area
//
/******************************************************************************/
void TriPanel::computeArea ()
{
  _area = tri_area(_verts[0]->x(), _verts[0]->y(), _verts[0]->z(),
                   _verts[1]->x(), _verts[1]->y(), _verts[1]->z(),
                   _verts[2]->x(), _verts[2]->y(), _verts[2]->z());
}

/******************************************************************************/
//
// Computes normal
//
/******************************************************************************/
void TriPanel::computeNormal ()
{
  _norm = tri_normal(_verts[0]->x(), _verts[0]->y(), _verts[0]->z(),
                     _verts[1]->x(), _verts[1]->y(), _verts[1]->z(),
                     _verts[2]->x(), _verts[2]->y(), _verts[2]->z());
}

/******************************************************************************/
//
// Computes centroid by intersection of vertex - edge midpoint lines
//
/******************************************************************************/
void TriPanel::computeCentroid () 
{
  _cen = tri_centroid(_verts[0]->x(), _verts[0]->y(), _verts[0]->z(),
                      _verts[1]->x(), _verts[1]->y(), _verts[1]->z(),
                      _verts[2]->x(), _verts[2]->y(), _verts[2]->z());
}

/******************************************************************************/
//
// Computes transform from inertial frame to panel frame and inverse transform.
// Also computes and stores panel endpoints in panel frame.
//
/******************************************************************************/
void TriPanel::computeTransform ()
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
    vec(0) = _verts[i]->x() - _cen[0];
    vec(1) = _verts[i]->y() - _cen[1];
    vec(2) = _verts[i]->z() - _cen[2];
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
double TriPanel::sourcePhiCoeff ( const double & x, const double & y, 
                                  const double & z, bool onpanel,
                                  const std::string & side,
                                  bool mirror_y ) const
{
  Eigen::Vector3d vec, transvec;
  double coeff;

  // Vector from centroid to point

  vec(0) = x - _cen[0];
  vec(1) = y - _cen[1];
  vec(2) = z - _cen[2];

  // Transform point to panel frame
  
  transvec = _trans*vec;
  
  // Farfield approximation

  if (vec.norm() > Panel::_farfield_distance_factor*_length)
    coeff = _area*point_source_potential(transvec(0), transvec(1),
                                         transvec(2));
  else
  {
    // Source potential -- panel endpoints given in clockwise order
    
    coeff = tri_source_potential(transvec(0), transvec(1), transvec(2),
                                 _xtrans[0], _ytrans[0], _xtrans[2], _ytrans[2],
                                 _xtrans[1], _ytrans[1], onpanel, side);
  }

  // Compute mirror image contribution if requested. Always assume the mirror
  // image point is not on the panel.

  if (mirror_y)
  {
    vec(0) =  x - _cen[0];
    vec(1) = -y - _cen[1];
    vec(2) =  z - _cen[2];
    transvec = _trans*vec;
    if (vec.norm() > Panel::_farfield_distance_factor*_length)
      coeff += _area*point_source_potential(transvec(0), transvec(1),
                                            transvec(2)); 
    else
      coeff += tri_source_potential(transvec(0), transvec(1), transvec(2),
                                 _xtrans[0], _ytrans[0], _xtrans[2], _ytrans[2],
                                 _xtrans[1], _ytrans[1], onpanel, side);
  }

  return coeff;
}

/******************************************************************************/
//
// Computes source induced velocity component
//
/******************************************************************************/
Eigen::Vector3d TriPanel::sourceVCoeff ( const double & x, const double & y, 
                                         const double & z, bool onpanel,
                                         const std::string & side,
                                         bool mirror_y ) const
{
  Eigen::Vector3d vec, transvec, velpf, velif, velif_mirror;

  // Vector from centroid to point

  vec(0) = x - _cen[0];
  vec(1) = y - _cen[1];
  vec(2) = z - _cen[2];

  // Transform point to panel frame

  transvec = _trans*vec;

  // Farfield approximation

  if (vec.norm() > Panel::_farfield_distance_factor*_length)
    velpf = _area*point_source_velocity(transvec(0), transvec(1), transvec(2));
  else
  {
    // Source velocity -- panel endpoints given in clockwise order
    
    velpf = tri_source_velocity(transvec(0), transvec(1), transvec(2),
                                _xtrans[0], _ytrans[0], _xtrans[2], _ytrans[2],
                                _xtrans[1], _ytrans[1], onpanel, side);
  }
  
  // Convert to inertial frame
  
  velif = _invtrans*velpf;

  // Compute mirror image contribution if requested. Always assume the mirror
  // image point is not on the panel.

  if (mirror_y)
  {
    vec(0) =  x - _cen[0];
    vec(1) = -y - _cen[1];
    vec(2) =  z - _cen[2];
    transvec = _trans*vec;
    if (vec.norm() > Panel::_farfield_distance_factor*_length)
      velpf = _area*point_source_velocity(transvec(0), transvec(1),
                                          transvec(2));
    else
      velpf = tri_source_velocity(transvec(0), transvec(1), transvec(2),
                                 _xtrans[0], _ytrans[0], _xtrans[2], _ytrans[2],
                                 _xtrans[1], _ytrans[1], false, side);
    velif_mirror = _invtrans*velpf;
    velif(0) += velif_mirror(0);
    velif(1) -= velif_mirror(1);
    velif(2) += velif_mirror(2);
  }

  return velif;
}

/******************************************************************************/
//
// Computes doublet influence coefficient
//
/******************************************************************************/
double TriPanel::doubletPhiCoeff ( const double & x, const double & y, 
                                   const double & z, bool onpanel,
                                   const std::string & side,
                                   bool mirror_y ) const
{
  Eigen::Vector3d vec, transvec;
  double coeff;

  // Transform point to panel frame

  vec(0) = x - _cen[0];
  vec(1) = y - _cen[1];
  vec(2) = z - _cen[2];
  transvec = _trans*vec;

  // Farfield approximation. Note that transformation is required here since
  // doublets are directional.

  if (vec.norm() > Panel::_farfield_distance_factor*_length)
    coeff = _area*point_doublet_potential(transvec(0), transvec(1),
                                          transvec(2));
  else
  {
    // Doublet potential -- points given in clockwise order
    
    coeff = tri_doublet_potential(transvec(0), transvec(1), transvec(2),
                                 _xtrans[0], _ytrans[0], _xtrans[2], _ytrans[2],
                                 _xtrans[1], _ytrans[1], onpanel, side);
  }

  // Compute mirror image contribution if requested. Always assume the mirror
  // image point is not on the panel.

  if (mirror_y)
  {
    vec(0) =  x - _cen[0];
    vec(1) = -y - _cen[1];
    vec(2) =  z - _cen[2];
    transvec = _trans*vec;
    if (vec.norm() > Panel::_farfield_distance_factor*_length)
      coeff += _area*point_doublet_potential(transvec(0), transvec(1),
                                             transvec(2));
    else
      coeff += tri_doublet_potential(transvec(0), transvec(1), transvec(2),
                                 _xtrans[0], _ytrans[0], _xtrans[2], _ytrans[2],
                                 _xtrans[1], _ytrans[1], onpanel, side);
  }

  return coeff;
}

/******************************************************************************/
//
// Computes doublet induced velocity component
//
/******************************************************************************/
Eigen::Vector3d TriPanel::doubletVCoeff ( const double & x, const double & y, 
                                         const double & z, bool onpanel,
                                         const std::string & side,
                                         bool mirror_y ) const
{
  Eigen::Vector3d vec, transvec, velpf, velif, velif_mirror;

  // Transform point to panel frame

  vec(0) = x - _cen[0];
  vec(1) = y - _cen[1];
  vec(2) = z - _cen[2];
  transvec = _trans*vec;

  // Farfield approximation. Note that transformation is required here since
  // doublets are directional.

  if (vec.norm() > Panel::_farfield_distance_factor*_length)
    velpf = _area*point_doublet_velocity(transvec(0), transvec(1), transvec(2));
  else
  {
    // Doublet velocity -- points given in clockwise order
    
    velpf = tri_doublet_velocity(transvec(0), transvec(1), transvec(2),
                                 _xtrans[0], _ytrans[0], _xtrans[2], _ytrans[2],
                                 _xtrans[1], _ytrans[1], onpanel, side);
  }

  // Convert to inertial frame

  velif = _invtrans*velpf;

  // Compute mirror image contribution if requested. Always assume the mirror
  // image point is not on the panel.

  if (mirror_y)
  {
    vec(0) =  x - _cen[0];
    vec(1) = -y - _cen[1];
    vec(2) =  z - _cen[2];
    transvec = _trans*vec;
    if (vec.norm() > Panel::_farfield_distance_factor*_length)
      velpf = _area*point_doublet_velocity(transvec(0), transvec(1),
                                           transvec(2));
    else
      velpf = tri_doublet_velocity(transvec(0), transvec(1), transvec(2),
                                 _xtrans[0], _ytrans[0], _xtrans[2], _ytrans[2],
                                 _xtrans[1], _ytrans[1], false, side);
    velif_mirror = _invtrans*velpf;
    velif(0) += velif_mirror(0);
    velif(1) -= velif_mirror(1);
    velif(2) += velif_mirror(2);
  }

  return velif;
}

/******************************************************************************/
//
// Computes induced velocity potential at a point
//
/******************************************************************************/
double TriPanel::inducedPotential ( const double & x, const double & y, 
                                    const double & z, bool onpanel,
                                    const std::string & side,
                                    bool mirror_y ) const
{
  double potential;

  potential = ( sourcePhiCoeff(x, y, z, onpanel, side, mirror_y)*_sigma
            +   doubletPhiCoeff(x, y, z, onpanel, side, mirror_y)*_mu );
  
  return potential;
}

/******************************************************************************/
//
// Computes induced velocity component at a point
//
/******************************************************************************/
Eigen::Vector3d TriPanel::inducedVelocity ( const double & x, const double & y, 
                                            const double & z, bool onpanel,
                                            const std::string & side,
                                            bool mirror_y ) const
{
  Eigen::Vector3d vel;

  vel = ( sourceVCoeff(x, y, z, onpanel, side, mirror_y)*_sigma
      +   doubletVCoeff(x, y, z, onpanel, side, mirror_y)*_mu );

  return vel;
}
