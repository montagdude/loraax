#include <vector>
#include <Eigen/Core>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "singularities.h"
#include "quadpanel.h"
#include "vertex.h"
#include "settings.h"

int main ()
{
  Vertex v1, v2, v3, v4;
  QuadPanel quadpanel;
  std::vector<std::vector<Vertex> > point_doublets;
  std::vector<double> y;
  std::vector<Eigen::Vector3d> velarr, velq;
  unsigned int npoints, i, j, k, nx, ny;
  double deltx, delty, pdoublet, px, py, x, z, ymin, ymax, dy;
  std::ofstream f;

  farfield_distance_factor = 1.E+06;

  v1.setCoordinates(0., -0.5, 0.0);
  v2.setCoordinates(0., 0.5, 0.0);
  v3.setCoordinates(1., 0.5, 0.);
  v4.setCoordinates(1., -0.5, 0.);

  quadpanel.addVertex(&v4);
  quadpanel.addVertex(&v3);
  quadpanel.addVertex(&v2);
  quadpanel.addVertex(&v1);
  quadpanel.setSourceStrength(0.0);
  quadpanel.setDoubletStrength(1.0);

  // Create an array of point doublets to compare with quad panel

  nx = 50;
  ny = 50;
  point_doublets.resize(nx);
  deltx = (1.0 - 0.0)/double(nx);
  delty = (0.5 - -0.5)/double(ny);
  pdoublet = 1.*(deltx*delty);		// Equivalent point doublet strength
  for ( j = 0; j < nx; j++ )
  {
    px = 0.0 + deltx/2. + double(j)*deltx;
    point_doublets[j].resize(ny);
    for ( k = 0; k < ny; k++ ) 
    {
      py = -0.5 + delty/2. + double(k)*delty;
      point_doublets[j][k].setCoordinates(px, py, 0.);
    }
  }

  npoints = 100;
  ymin = -5;
  ymax = 5;
  dy = (ymax-ymin)/double(npoints-1);
  y.resize(npoints);
  velarr.resize(npoints);
  velq.resize(npoints);
  x = 0.5;
  z = 1.0;
  for ( i = 0; i < npoints; i++ )
  {
    y[i] = ymin + double(i)*dy;
    velarr[i](0) = 0.;
    velarr[i](1) = 0.;
    velarr[i](2) = 0.;
    for ( j = 0; j < nx; j++ )
    {
      for ( k = 0; k < ny; k++ )
      {
        velarr[i] += point_doublet_velocity(x-point_doublets[j][k].x(),
                                           y[i]-point_doublets[j][k].y(),
                                           z-point_doublets[j][k].z())*pdoublet;
      }
    }
    velq[i] = quadpanel.inducedVelocity(x, y[i], z, false, "top");
  } 

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

  f.open("doubletarr1.dat");
  f << "VARIABLES = \"x\" \"y\" \"z\" \"Vx\" \"Vy\" \"Vz\"" << std::endl;
  for ( i = 0; i < npoints; i++ )
  {
    f << std::setprecision(9) << x << " " << y[i] << " " << z << " "
      << velarr[i](0) << " " << velarr[i](1) << " " << velarr[i](2) << " "
      << std::endl;
  }
  f.close();
  std::cout << "Wrote file doubletarr1.dat." << std::endl;

  return 0;
}
