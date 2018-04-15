#include <Eigen/Core>
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
  Eigen::Vector3d velv, velq;
  double x, y, z, rcore;

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

  x = 0.5;
  y = 0.2;
  z = 0.0;
  velv = vring.inducedVelocity(x, y, z, rcore);
  velq = quadpanel.inducedVelocity(x, y, z, true, "top");

  std::cout << "Vortex velocity: " << std::setprecision(9) << velv(0) << " "
            << velv(1) << " " << velv(2) << std::endl;
  std::cout << "Quad velocity: " << std::setprecision(9) << velq(0) << " "
            << velq(1) << " " << velq(2) << std::endl;

  return 0;
}
