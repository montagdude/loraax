#include <vector>
#include <Eigen/Core>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "vortex_ring.h"
#include "vertex.h"
#include "settings.h"

int main ()
{
  Vertex v1, v2, v3, v4;
  VortexRing vring;
  std::vector<double> y, vmag;  
  std::vector<Eigen::Vector3d> vel;
  unsigned int npoints, i;
  double x, z, ymin, ymax, dy, rcore;
  std::ofstream f;

  farfield_distance_factor = 1.E+06;

  v1.setCoordinates(0., -0.5, 0.);
  v2.setCoordinates(0., 0.5, 0.);
  v3.setCoordinates(1., 0.5, 0.);
  v4.setCoordinates(1., -0.5, 0.);
  vring.addVertex(&v1);
  vring.addVertex(&v2);
  vring.addVertex(&v3);
  vring.addVertex(&v4);
  vring.setCirculation(1.0);
 
  npoints = 100;
  ymin = -5;
  ymax = 5;
  dy = (ymax-ymin)/double(npoints-1);
  y.resize(npoints);
  vel.resize(npoints);
  vmag.resize(npoints);
  rcore = 1.E-08;
  x = 0.5;
  z = 1.0;
  for ( i = 0; i < npoints; i++ )
  {
    y[i] = ymin + double(i)*dy;
    vel[i] = vring.inducedVelocity(x, y[i], z, rcore);
    vmag[i] = vel[i].norm();
  } 

  f.open("vortex_ring_data.dat");
  f << "VARIABLES = \"x\" \"y\" \"z\" \"Vx\" \"Vy\" \"Vz\" \"Vmag\"" << std::endl;
  for ( i = 0; i < npoints; i++ )
  {
    f << std::setprecision(9) << x << " " << y[i] << " " << z << " "
      << vel[i](0) << " " << vel[i](1) << " " << vel[i](2) << " " << vmag[i]
      << std::endl;
  }
  f.close();
  std::cout << "Wrote file vortex_ring_data.dat." << std::endl;

  return 0;
}
