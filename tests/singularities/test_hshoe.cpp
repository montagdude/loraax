#include <vector>
#include <Eigen/Core>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "vortex_ring.h"
#include "horseshoe_vortex.h"
#include "quadpanel.h"
#include "vertex.h"
#include "settings.h"

int main ()
{
  Vertex v1, v2, v3, v4;
  VortexRing vring;
  HorseshoeVortex hshoe;
  QuadPanel pan;
  std::vector<double> y;
  std::vector<Eigen::Vector3d> velv, velq, velh;
  unsigned int npoints, i;
  double x, z, ymin, ymax, dy, rcore;
  std::ofstream f;

  v1.setCoordinates(0., -0.5, 0.);
  v2.setCoordinates(0., 0.5, 0.);
  v3.setCoordinates(1., 0.5, 0.);
  v4.setCoordinates(1., -0.5, 0.);

  vring.addVertex(&v4);
  vring.addVertex(&v1);
  vring.addVertex(&v2);
  vring.addVertex(&v3);
  vring.setCirculation(1.0);

  pan.addVertex(&v4);
  pan.addVertex(&v3);
  pan.addVertex(&v2);
  pan.addVertex(&v1);
  pan.setSourceStrength(0.0);
  pan.setDoubletStrength(1.0);
 
  hshoe.addVertex(&v4);
  hshoe.addVertex(&v1);
  hshoe.addVertex(&v2);
  hshoe.addVertex(&v3);
  hshoe.setCirculation(1.0);

  npoints = 100;
  ymin = -5.;
  ymax = 5;
  dy = (ymax-ymin)/double(npoints-1);
  y.resize(npoints);
  velv.resize(npoints);
  velq.resize(npoints);
  velh.resize(npoints);
  rcore = 1.E-12;
  x = 0.5;
  z = 1.0;
  for ( i = 0; i < npoints; i++ )
  {
    y[i] = ymin + double(i)*dy;
    velv[i] = vring.inducedVelocity(x, y[i], z, rcore);
    velq[i] = pan.inducedVelocity(x, y[i], z, false, "top", false);
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

  f.open("quad_doublet3.dat");
  f << "VARIABLES = \"x\" \"y\" \"z\" \"Vx\" \"Vy\" \"Vz\"" << std::endl;
  for ( i = 0; i < npoints; i++ )
  {
    f << std::setprecision(9) << x << " " << y[i] << " " << z << " "
      << velq[i](0) << " " << velq[i](1) << " " << velq[i](2) << " "
      << std::endl;
  }
  f.close();
  std::cout << "Wrote file quad_doublet3.dat." << std::endl;

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

  // Move last two vertices back and recalculate vortex ring and quad panel
  // (horseshoe vortex result stays the same as it should)

  v3.translate(1., 0., 0.);
  v4.translate(1., 0., 0.);
  pan.recomputeGeometry();
  for ( i = 0; i < npoints; i++ )
  {
    y[i] = ymin + double(i)*dy;
    velv[i] = vring.inducedVelocity(x, y[i], z, rcore);
    velq[i] = pan.inducedVelocity(x, y[i], z, false, "top", false);
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

  f.open("quad_doublet4.dat");
  f << "VARIABLES = \"x\" \"y\" \"z\" \"Vx\" \"Vy\" \"Vz\"" << std::endl;
  for ( i = 0; i < npoints; i++ )
  {
    f << std::setprecision(9) << x << " " << y[i] << " " << z << " "
      << velq[i](0) << " " << velq[i](1) << " " << velq[i](2) << " "
      << std::endl;
  }
  f.close();
  std::cout << "Wrote file quad_doublet4.dat." << std::endl;

  // Do it one more time

  v3.translate(8., 0., 0.);
  v4.translate(8., 0., 0.);
  pan.recomputeGeometry();
  for ( i = 0; i < npoints; i++ )
  {
    y[i] = ymin + double(i)*dy;
    velv[i] = vring.inducedVelocity(x, y[i], z, rcore);
    velq[i] = pan.inducedVelocity(x, y[i], z, false, "top", false);
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

  f.open("quad_doublet5.dat");
  f << "VARIABLES = \"x\" \"y\" \"z\" \"Vx\" \"Vy\" \"Vz\"" << std::endl;
  for ( i = 0; i < npoints; i++ )
  {
    f << std::setprecision(9) << x << " " << y[i] << " " << z << " "
      << velq[i](0) << " " << velq[i](1) << " " << velq[i](2) << " "
      << std::endl;
  }
  f.close();
  std::cout << "Wrote file quad_doublet5.dat." << std::endl;

  return 0;
}
