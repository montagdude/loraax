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
  QuadPanel pan1;
  std::vector<double> y;
  std::vector<Eigen::Vector3d> velv, velq;
  unsigned int npoints, i;
  double x, z, ymin, ymax, dy, rcore;
  std::ofstream f;

  v1.setCoordinates(0., 0.0, 0.0);
  v2.setCoordinates(1., 0.0, 0.0);
  v3.setCoordinates(1., 0.5, 0.5);
  v4.setCoordinates(0., 0.5, 0.5);

  vring.addVertex(&v4);
  vring.addVertex(&v3);
  vring.addVertex(&v2);
  vring.addVertex(&v1);
  vring.setCirculation(1.0);
 
  pan1.addVertex(&v4);
  pan1.addVertex(&v1);
  pan1.addVertex(&v2);
  pan1.addVertex(&v3);
  pan1.setSourceStrength(0.0);
  pan1.setDoubletStrength(1.0);

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
    velq[i] = pan1.inducedVelocity(x, y[i], z, false, "top");
  } 

  f.open("vring5.dat");
  f << "VARIABLES = \"x\" \"y\" \"z\" \"Vx\" \"Vy\" \"Vz\"" << std::endl;
  for ( i = 0; i < npoints; i++ )
  {
    f << std::setprecision(9) << x << " " << y[i] << " " << z << " "
      << velv[i](0) << " " << velv[i](1) << " " << velv[i](2) << " "
      << std::endl;
  }
  f.close();
  std::cout << "Wrote file vring5.dat." << std::endl;

  f.open("quad_doublet1.dat");
  f << "VARIABLES = \"x\" \"y\" \"z\" \"Vx\" \"Vy\" \"Vz\"" << std::endl;
  for ( i = 0; i < npoints; i++ )
  {
    f << std::setprecision(9) << x << " " << y[i] << " " << z << " "
      << velq[i](0) << " " << velq[i](1) << " " << velq[i](2) << " "
      << std::endl;
  }
  f.close();
  std::cout << "Wrote file quad_doublet1.dat." << std::endl;

  return 0;
}
