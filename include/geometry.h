// Contains functions to compute geometric quantities (e.g. area, normal,
// centroid)

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <Eigen/Core>

// Public functions

double tri_area ( const double &, const double &, const double &,
                  const double &, const double &, const double &,
                  const double &, const double &, const double & );
							// Computes area of triangle from three nodes
							//   given in counterclockwise order

double quad_area ( const double &, const double &, const double &,
                   const double &, const double &, const double &,
                   const double &, const double &, const double &,
                   const double &, const double &, const double & );
							// Computes area of quadrilateral from four nodes
							//   given in counterclockwise order

Eigen::Vector3d tri_normal ( const double &, const double &, const double &,
                             const double &, const double &, const double &,
                             const double &, const double &, const double & );
							// Computes normal of triangle from three nodes
							//   given in counterclockwise order

Eigen::Vector3d quad_normal ( const double &, const double &, const double &,
                              const double &, const double &, const double &,
                              const double &, const double &, const double &,
                              const double &, const double &, const double & );
							// Computes normal of quadrilateral from four nodes
							//   given in counterclockwise order

Eigen::Vector3d tri_centroid ( const double &, const double &, const double &,
                               const double &, const double &, const double &,
                               const double &, const double &, const double & );
							// Computes centroid of triangle from three nodes
							//   given in counterclockwise order via intersection
							//   of midpoint lines

Eigen::Vector3d quad_centroid ( const double &, const double &, const double &,
                                const double &, const double &, const double &,
                                const double &, const double &, const double &,
                                const double &, const double &, const double & 
                              );
							// Computes centroid of quadrilateral from four
							//   nodes given in counterclockwise order by
							//   breaking up into two triangles

Eigen::Matrix3d transform_from_normal ( const double &, const double &,
                                        const double & );
							// Computes inertial->panel transformation from
							//   normal vector

std::vector<double> uniform_spacing ( const double & slen, unsigned int n );
std::vector<double> cosine_spacing ( const double & slen, unsigned int n );
std::vector<double> sine_spacing ( const double & slen, unsigned int n );
							// Edge spacing functions

void compute_plane ( const Eigen::Vector3d & point,
                     const Eigen::Vector3d & norm, double & a, double & b,
                     double & c, double & d );
							// Equation of a plane from point and normal

Eigen::Vector3d line_plane_intersection ( const Eigen::Vector3d & p0,
                                    const Eigen::Vector3d & l, const double & a,
                                    const double & b, const double & c,
                                    const double & d );
#endif
