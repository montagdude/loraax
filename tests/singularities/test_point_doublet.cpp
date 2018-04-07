#include <Eigen/Core>
#include <iostream>
#include <iomanip>
#include "singularities.h"

int main ()
{
  double x, y, z;
  double h;
  Eigen::Vector3d vel, velapprox;
  double mu;

  x = 0.1;
  y = -0.1;
  z = 0.2;
  mu = 1.0; 

  // Velocity from analytical expression

  vel = point_doublet_velocity(x, y, z)*mu;

  // Velocity from finite differences

  h = 1.e-10;
  velapprox(0) = ( point_doublet_potential(x+h/2., y, z)
               -   point_doublet_potential(x-h/2., y, z) ) / h;
  velapprox(1) = ( point_doublet_potential(x, y+h/2., z)
               -   point_doublet_potential(x, y-h/2., z) ) / h;
  velapprox(2) = ( point_doublet_potential(x, y, z+h/2.)
               -   point_doublet_potential(x, y, z-h/2.) ) / h;

  std::cout << "Analytical    Finite differences" << std::endl;
  std::cout << std::setprecision(9) << vel(0) << "  " << velapprox(0)
            << std::endl;
  std::cout << std::setprecision(9) << vel(1) << "  " << velapprox(1)
            << std::endl;
  std::cout << std::setprecision(9) << vel(2) << "  " << velapprox(2)
            << std::endl;

  return 0;
}
