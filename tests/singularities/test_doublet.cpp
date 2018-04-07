#include <vector>
#include <Eigen/Core>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "vortex_ring.h"
#include "quadpanel.h"
#include "vertex.h"
#include "settings.h"

int main ()
{
  Vertex v1, v2, v3, v4;
  VortexRing vring;
  QuadPanel quadpanel;
  std::vector<double> y;
  std::vector<Eigen::Vector3d> velv, velq;
  unsigned int npoints, i;
  double x, z, ymin, ymax, dy, rcore;
  std::ofstream f;

  farfield_distance_factor = 1.E+06;

  v1.setCoordinates(0., -0.5, 0.0);
  v2.setCoordinates(0., 0.5, 0.0);
  v3.setCoordinates(1., 0.5, 0.);
  v4.setCoordinates(1., -0.5, 0.);

  vring.addVertex(&v4);
  vring.addVertex(&v1);
  vring.addVertex(&v2);
  vring.addVertex(&v3);
  vring.setCirculation(1.0);
 
  quadpanel.addVertex(&v4);
  quadpanel.addVertex(&v3);
  quadpanel.addVertex(&v2);
  quadpanel.addVertex(&v1);
  quadpanel.setSourceStrength(0.0);
  quadpanel.setDoubletStrength(1.0);

  npoints = 100;
  ymin = -5;
  ymax = 5;
  dy = (ymax-ymin)/double(npoints-1);
  y.resize(npoints);
  velv.resize(npoints);
  velq.resize(npoints);
  rcore = 1.E-12;
  x = 0.5;
  z = 1.0;
  for ( i = 0; i < npoints; i++ )
  {
    y[i] = ymin + double(i)*dy;
    velv[i] = vring.inducedVelocity(x, y[i], z, rcore);
    velq[i] = quadpanel.inducedVelocity(x, y[i], z, false, "top");
  } 

  f.open("vring1.dat");
  f << "VARIABLES = \"x\" \"y\" \"z\" \"Vx\" \"Vy\" \"Vz\"" << std::endl;
  for ( i = 0; i < npoints; i++ )
  {
    f << std::setprecision(9) << x << " " << y[i] << " " << z << " "
      << velv[i](0) << " " << velv[i](1) << " " << velv[i](2) << " "
      << std::endl;
  }
  f.close();
  std::cout << "Wrote file vring1.dat." << std::endl;

  f.open("doublet1.dat");
  f << "VARIABLES = \"x\" \"y\" \"z\" \"Vx\" \"Vy\" \"Vz\"" << std::endl;
  for ( i = 0; i < npoints; i++ )
  {
    f << std::setprecision(9) << x << " " << y[i] << " " << z << " "
      << velq[i](0) << " " << velq[i](1) << " " << velq[i](2) << " "
      << std::endl;
  }
  f.close();
  std::cout << "Wrote file doublet1.dat." << std::endl;

  return 0;
}
