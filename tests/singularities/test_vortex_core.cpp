#include <vector>
#include <Eigen/Core>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "singularities.h"

int main ()
{
  std::vector<double> y;
  std::vector<Eigen::Vector3d> velv1, velv2;
  unsigned int npoints, i;
  double x1, x2, y1, y2, z1, z2;
  double x, z, ymin, ymax, dy;
  const double eps = 1e-12;
  std::ofstream f;

  x1 = 0.;
  x2 = 5.;
  y1 = 0.;
  y2 = 0.;
  z1 = 0.;
  z2 = 0.;

  npoints = 101;
  ymin = -0.1;
  ymax = 0.1;
  dy = (ymax-ymin)/double(npoints-1);
  y.resize(npoints);
  velv1.resize(npoints);
  velv2.resize(npoints);
  x = 2.5;
  z = 0.0;
  for ( i = 0; i < npoints; i++ )
  {
    y[i] = ymin + double(i)*dy;
    velv1[i] = vortex_velocity(x, y[i], z, x1, y1, z1, x2, y2, z2, 1e-3, false);
    velv2[i] = vortex_velocity(x, y[i], z, x1, y1, z1, x2, y2, z2, eps, false);
  } 

  f.open("vortex_core1.dat");
  f << "VARIABLES = \"x\" \"y\" \"z\" \"Vx\" \"Vy\" \"Vz\"" << std::endl;
  for ( i = 0; i < npoints; i++ )
  {
    f << std::setprecision(9) << x << " " << y[i] << " " << z << " "
      << velv1[i](0) << " " << velv1[i](1) << " " << velv1[i](2) << " "
      << std::endl;
  }
  f.close();
  std::cout << "Wrote file vortex_core1.dat." << std::endl;

  f.open("vortex_core2.dat");
  f << "VARIABLES = \"x\" \"y\" \"z\" \"Vx\" \"Vy\" \"Vz\"" << std::endl;
  for ( i = 0; i < npoints; i++ )
  {
    f << std::setprecision(9) << x << " " << y[i] << " " << z << " "
      << velv2[i](0) << " " << velv2[i](1) << " " << velv2[i](2) << " "
      << std::endl;
  }
  f.close();
  std::cout << "Wrote file vortex_core2.dat." << std::endl;

  return 0;
}
