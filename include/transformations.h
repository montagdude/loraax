// Contains functions to perform linear transformations

#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H

#include <Eigen/Core>
#include <string>

// Public routines

Eigen::Matrix3d euler_rotation (
                       const double & phid, const double & thetad, 
                       const double & psid, const std::string & order = "321" );
Eigen::Matrix3d inverse_euler_rotation (
                       const double & phid, const double & thetad, 
                       const double & psid, const std::string & order = "321" );

#endif
