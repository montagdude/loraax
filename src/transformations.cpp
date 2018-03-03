// Contains functions to perform linear transformations

#include <cmath>
#include <Eigen/Core>

/******************************************************************************/
//
// 3-2-1 Euler rotation
//
/******************************************************************************/
Eigen::Matrix3d euler_rotation(const double & phi, const double & theta, 
                               const double & psi)
{
  Eigen::Matrix3d L1, L2, L3, transform;

  L1.setIdentity();
  L2.setIdentity();
  L3.setIdentity();

  L1(1,1) = cos(phi);  L1(1,2) = sin(phi);
  L1(2,1) = -sin(phi); L1(2,2) = cos(phi);

  L2(0,0) = cos(theta); L2(0,2) = -sin(theta);
  L2(2,0) = sin(theta); L2(2,2) = cos(theta);

  L3(0,0) = cos(psi);  L3(0,1) = sin(psi);
  L3(1,0) = -sin(psi); L3(1,1) = cos(psi);

  transform = L1*L2*L3;
  return transform;
}
