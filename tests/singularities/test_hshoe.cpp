#include <vector>
#include <Eigen/Core>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "vortex_ring.h"
#include "horseshoe_vortex.h"
#include "vertex.h"
#include "settings.h"

int main ()
{
  Vertex v1, v2, v3, v4;
  VortexRing vring;
  HorseshoeVortex hshoe;
  std::vector<double> y;
  std::vector<Eigen::Vector3d> velv, velh;
  unsigned int npoints, i;
  double x, z, ymin, ymax, dy, rcore;
  std::ofstream f;

  farfield_distance_factor = 1.E+06;

  v1.setCoordinates(0., -0.5, 0.);
  v2.setCoordinates(0., 0.5, 0.);
  v3.setCoordinates(1., 0.5, 0.);
  v4.setCoordinates(1., -0.5, 0.);

  vring.addVertex(&v4);
  vring.addVertex(&v1);
  vring.addVertex(&v2);
  vring.addVertex(&v3);
  vring.setCirculation(1.0);
 
  hshoe.addVertex(&v4);
  hshoe.addVertex(&v1);
  hshoe.addVertex(&v2);
  hshoe.addVertex(&v3);
  hshoe.setCirculation(1.0);

  npoints = 100;
  ymin = -5;
  ymax = 5;
  dy = (ymax-ymin)/double(npoints-1);
  y.resize(npoints);
  velv.resize(npoints);
  velh.resize(npoints);
  rcore = 1.E-12;
  x = 0.5;
  z = 1.0;
  for ( i = 0; i < npoints; i++ )
  {
    y[i] = ymin + double(i)*dy;
    velv[i] = vring.inducedVelocity(x, y[i], z, rcore);
    velh[i] = hshoe.inducedVelocity(x, y[i], z, rcore);
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

  f.open("hshoe1.dat");
  f << "VARIABLES = \"x\" \"y\" \"z\" \"Vx\" \"Vy\" \"Vz\"" << std::endl;
  for ( i = 0; i < npoints; i++ )
  {
    f << std::setprecision(9) << x << " " << y[i] << " " << z << " "
      << velh[i](0) << " " << velh[i](1) << " " << velh[i](2) << " "
      << std::endl;
  }
  f.close();
  std::cout << "Wrote file hshoe1.dat." << std::endl;

  // Move last two vertices back and recalculate vortex ring
  // (horseshoe vortex result stays the same as it should)

  v3.translate(1., 0., 0.);
  v4.translate(1., 0., 0.);
  for ( i = 0; i < npoints; i++ )
  {
    y[i] = ymin + double(i)*dy;
    velv[i] = vring.inducedVelocity(x, y[i], z, rcore);
    velh[i] = hshoe.inducedVelocity(x, y[i], z, rcore);
  } 

  f.open("vring2.dat");
  f << "VARIABLES = \"x\" \"y\" \"z\" \"Vx\" \"Vy\" \"Vz\"" << std::endl;
  for ( i = 0; i < npoints; i++ )
  {
    f << std::setprecision(9) << x << " " << y[i] << " " << z << " "
      << velv[i](0) << " " << velv[i](1) << " " << velv[i](2) << " "
      << std::endl;
  }
  f.close();
  std::cout << "Wrote file vring2.dat." << std::endl;

  // Do it one more time

  v3.translate(8., 0., 0.);
  v4.translate(8., 0., 0.);
  for ( i = 0; i < npoints; i++ )
  {
    y[i] = ymin + double(i)*dy;
    velv[i] = vring.inducedVelocity(x, y[i], z, rcore);
    velh[i] = hshoe.inducedVelocity(x, y[i], z, rcore);
  } 

  f.open("vring3.dat");
  f << "VARIABLES = \"x\" \"y\" \"z\" \"Vx\" \"Vy\" \"Vz\"" << std::endl;
  for ( i = 0; i < npoints; i++ )
  {
    f << std::setprecision(9) << x << " " << y[i] << " " << z << " "
      << velv[i](0) << " " << velv[i](1) << " " << velv[i](2) << " "
      << std::endl;
  }
  f.close();
  std::cout << "Wrote file vring3.dat." << std::endl;

  return 0;
}
