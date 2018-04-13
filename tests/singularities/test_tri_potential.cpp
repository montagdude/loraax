#include <vector>
#include <Eigen/Core>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "tripanel.h"
#include "vertex.h"
#include "settings.h"

int main ()
{
  Vertex v1, v2, v3;
  TriPanel tripanel;
  std::vector<double> y, yfd;
  std::vector<Eigen::Vector3d> velexact, velfd;
  unsigned int npoints, i;
  double x, z, ymin, ymax, dy, phip, phim;
  std::ofstream f;

  v1.setCoordinates(0., 0.0, -0.2);
  v2.setCoordinates(1., 0.3, -0.2);
  v3.setCoordinates(0.5, 0.5, 0.5);

  tripanel.addVertex(&v3);
  tripanel.addVertex(&v2);
  tripanel.addVertex(&v1);
  tripanel.setSourceStrength(1.0);
  tripanel.setDoubletStrength(1.0);

  npoints = 100;
  ymin = -1;
  ymax = 1;
  dy = (ymax-ymin)/double(npoints-1);
  y.resize(npoints);
  velexact.resize(npoints);
  yfd.resize(npoints-1);
  velfd.resize(npoints-1);
  x = 0.3;
  z = 1.0;

  // This test has both a source and doublet and includes a mirror image panel.
  // The purpose is to test that the exact induced velocity equations match the
  // numerical gradient of the potential.

  for ( i = 0; i < npoints; i++ )
  {
    y[i] = ymin + double(i)*dy;
    velexact[i] = tripanel.inducedVelocity(x, y[i], z, false, "top", true);
  } 

  for ( i = 1; i < npoints; i++ )
  {
    yfd[i-1] = 0.5*(y[i-1] + y[i]);
    phip = tripanel.inducedPotential(x+dy/2., yfd[i-1], z, false, "top", true);
    phim = tripanel.inducedPotential(x-dy/2., yfd[i-1], z, false, "top", true);
    velfd[i-1](0) = (phip-phim) / dy;
    phip = tripanel.inducedPotential(x, y[i], z, false, "top", true);
    phim = tripanel.inducedPotential(x, y[i-1], z, false, "top", true);
    velfd[i-1](1) = (phip-phim) / dy;
    phip = tripanel.inducedPotential(x, yfd[i-1], z+dy/2., false, "top", true);
    phim = tripanel.inducedPotential(x, yfd[i-1], z-dy/2., false, "top", true);
    velfd[i-1](2) = (phip-phim) / dy;
  }

  f.open("tri_panel1.dat");
  f << "VARIABLES = \"x\" \"y\" \"z\" \"Vx\" \"Vy\" \"Vz\"" << std::endl;
  for ( i = 0; i < npoints; i++ )
  {
    f << std::setprecision(9) << x << " " << y[i] << " " << z << " "
      << velexact[i](0) << " " << velexact[i](1) << " " << velexact[i](2) << " "
      << std::endl;
  }
  f.close();
  std::cout << "Wrote file tri_panel1.dat." << std::endl;

  f.open("tri_panel2.dat");
  f << "VARIABLES = \"x\" \"y\" \"z\" \"Vx\" \"Vy\" \"Vz\"" << std::endl;
  for ( i = 0; i < npoints-1; i++ )
  {
    f << std::setprecision(9) << x << " " << yfd[i] << " " << z << " "
      << velfd[i](0) << " " << velfd[i](1) << " " << velfd[i](2) << " "
      << std::endl;
  }
  f.close();
  std::cout << "Wrote file tri_panel2.dat." << std::endl;

  return 0;
}
