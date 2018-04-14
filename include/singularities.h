// Computes velocity potential and velocity due to singularity elements

#ifndef SINGULARITIES_H
#define SINGULARITIES_H

#include <string>
#include <Eigen/Core>

// Public routines

// Note: all distributed source/doublet distributions on panels assume that 
// panel centroid is at the origin and that the panel endpoints lie on the z = 0
// plane.

double point_source_potential ( const double &, const double &, 
                                const double & );
                        // Source potential at a point due to point source at
                        //   the origin with unit strength
                       
double tri_source_potential ( const double &, const double &, const double &,
                              const double &, const double &, const double &, 
                              const double &, const double &, const double &, 
                              bool, const std::string & );
                        // Source potential at a point due to triangular
                        //   source panel at z = 0 with unit strength

double quad_source_potential ( const double &, const double &, const double &,
                               const double &, const double &, const double &, 
                               const double &, const double &, const double &, 
                               const double &, const double &, bool,
                               const std::string & );
                        // Source potential at a point due to quadrilateral
                        //   source panel at z = 0 with unit strength
                          
Eigen::Vector3d point_source_velocity ( const double &, const double &, 
                                        const double & );
                        // Source velocity at a point due to point source at
                        //   the origin with unit strength

Eigen::Vector3d tri_source_velocity ( const double &, const double &, 
                                      const double &, const double &, 
                                      const double &, const double &, 
                                      const double &, const double &, 
                                      const double &, bool, 
                                      const std::string & );
                        // Source velocity at a point due to triangular
                        //   source panel at z = 0 with unit strength

Eigen::Vector3d quad_source_velocity ( const double &, const double &, 
                                       const double &, const double &, 
                                       const double &, const double &, 
                                       const double &, const double &, 
                                       const double &, const double &, 
                                       const double &, bool, 
                                       const std::string & );
                        // Source velocity at a point due to quadrilateral
                        //   source panel at z = 0 with unit strength

double point_doublet_potential ( const double &, const double &, 
                                 const double & );
                        // Doublet potential at a point due to point doublet at
                        //   the origin with unit strength
                       
double tri_doublet_potential ( const double &, const double &, const double &,
                               const double &, const double &, const double &, 
                               const double &, const double &, const double &, 
                               bool, const std::string & );
                        // Doublet potential at a point due to triangular
                        //   doublet panel at z = 0 with unit strength

double quad_doublet_potential ( const double &, const double &, const double &,
                                const double &, const double &, const double &, 
                                const double &, const double &, const double &, 
                                const double &, const double &, bool,
                                const std::string & );
                        // Doublet potential at a point due to quadrilateral
                        //   doublet panel at z = 0 with unit strength

Eigen::Vector3d point_doublet_velocity ( const double &, const double &, 
                                         const double & );
                        // Doublet velocity at a point due to point doublet at
                        //   the origin with unit strength

Eigen::Vector3d tri_doublet_velocity ( const double &, const double &, 
                                       const double &, const double &, 
                                       const double &, const double &, 
                                       const double &, const double &, 
                                       const double &, bool, 
                                       const std::string & );
                        // Doublet velocity at a point due to triangular
                        //   doublet panel at z = 0 with unit strength

Eigen::Vector3d quad_doublet_velocity ( const double &, const double &, 
                                        const double &, const double &, 
                                        const double &, const double &, 
                                        const double &, const double &, 
                                        const double &, const double &, 
                                        const double &, bool, 
                                        const std::string & );
                        // Doublet velocity at a point due to quadrilateral
                        //   doublet panel at z = 0 with unit strength

Eigen::Vector3d vortex_velocity ( const double &, const double &, 
                                  const double &, const double &, 
                                  const double &, const double &,
                                  const double &, const double &, 
                                  const double &, const double &, bool );
                        // Velocity at a point due to finite- or semi-infinite-
                        //   length vortex filament with finite core and unit 
                        //   strength
#endif
