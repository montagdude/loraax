// Contains functions to perform linear transformations

#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H

#include <Eigen/Core>

// Public routines

Eigen::Matrix3d euler_rotation ( const double &, const double &, 
                                 const double & );

#endif
