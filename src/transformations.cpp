// Contains functions to perform linear transformations

#define _USE_MATH_DEFINES

#include <cmath>
#include <string>
#include <Eigen/Core>

/******************************************************************************/
//
// Euler rotation (inertial->body), angles given in degrees
//
/******************************************************************************/
Eigen::Matrix3d euler_rotation (const double & phid, const double & thetad, 
                                const double & psid, const std::string & order)
{
  Eigen::Matrix3d L1, L2, L3, transform;
  double phi, theta, psi, d2r;

  d2r = M_PI/180.;
  phi = phid*d2r;
  theta = thetad*d2r;
  psi = psid*d2r;

  L1.setIdentity();
  L2.setIdentity();
  L3.setIdentity();

  L1(1,1) = cos(phi);  L1(1,2) = sin(phi);
  L1(2,1) = -sin(phi); L1(2,2) = cos(phi);

  L2(0,0) = cos(theta); L2(0,2) = -sin(theta);
  L2(2,0) = sin(theta); L2(2,2) = cos(theta);

  L3(0,0) = cos(psi);  L3(0,1) = sin(psi);
  L3(1,0) = -sin(psi); L3(1,1) = cos(psi);

  if (order == "123")
    transform = L3*L2*L1;
  else
    transform = L1*L2*L3;

  return transform;
}

/******************************************************************************/
//
// Inverse euler rotation (body->inertial)
//
/******************************************************************************/
Eigen::Matrix3d inverse_euler_rotation (
                                const double & phid, const double & thetad,
                                const double & psid, const std::string & order )
{
  return euler_rotation(phid, thetad, psid, order).transpose();
}
