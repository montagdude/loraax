#include <vector>
#include <Eigen/Core>
#include <fstream>
#include <iostream>
#include <iomanip>
#include "singularities.h"
#include "quadpanel.h"
#include "tripanel.h"
#include "vertex.h"
#include "settings.h"

int main ()
{
  Vertex v1, v2, v3, v4;
  QuadPanel quadpanel;
  TriPanel tri1, tri2;
  std::vector<std::vector<Vertex> > point_sources;
  std::vector<double> x;
  std::vector<Eigen::Vector3d> velarr, velq, velt;
  unsigned int npoints, i, j, k, nx, ny;
  double deltx, delty, psource, px, py, y, z, xmin, xmax, dx;
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
  quadpanel.setSourceStrength(1.0);
  quadpanel.setDoubletStrength(0.0);

  tri1.addVertex(&v4);
  tri1.addVertex(&v3);
  tri1.addVertex(&v2);
  tri1.setSourceStrength(1.0);
  tri1.setDoubletStrength(0.0);

  tri2.addVertex(&v2);
  tri2.addVertex(&v1);
  tri2.addVertex(&v4);
  tri2.setSourceStrength(1.0);
  tri2.setDoubletStrength(0.0);

  // Create an array of point sources to compare with quad panel

  nx = 50;
  ny = 50;
  point_sources.resize(nx);
  deltx = (1.0 - 0.0)/double(nx);
  delty = (0.5 - -0.5)/double(ny);
  psource = 1.*(deltx*delty);		// Equivalent point source strength
  for ( j = 0; j < nx; j++ )
  {
    px = 0.0 + deltx/2. + double(j)*deltx;
    point_sources[j].resize(ny);
    for ( k = 0; k < ny; k++ ) 
    {
      py = -0.5 + delty/2. + double(k)*delty;
      point_sources[j][k].setCoordinates(px, py, 0.);
    }
  }

  npoints = 100;
  xmin = -5;
  xmax = 5;
  dx = (xmax-xmin)/double(npoints-1);
  x.resize(npoints);
  velarr.resize(npoints);
  velq.resize(npoints);
  velt.resize(npoints);
  y = 0.5;
  z = -1.0;
  for ( i = 0; i < npoints; i++ )
  {
    x[i] = xmin + double(i)*dx;
    velarr[i](0) = 0.;
    velarr[i](1) = 0.;
    velarr[i](2) = 0.;
    for ( j = 0; j < nx; j++ )
    {
      for ( k = 0; k < ny; k++ )
      {
        velarr[i] += point_source_velocity(x[i]-point_sources[j][k].x(),
                                           y-point_sources[j][k].y(),
                                           z-point_sources[j][k].z())*psource;
      }
    }
    velq[i] = quadpanel.inducedVelocity(x[i], y, z, false, "top");
    velt[i] = tri1.inducedVelocity(x[i], y, z, false, "top")
            + tri2.inducedVelocity(x[i], y, z, false, "top");
  } 

  f.open("source1.dat");
  f << "VARIABLES = \"x\" \"y\" \"z\" \"Vx\" \"Vy\" \"Vz\"" << std::endl;
  for ( i = 0; i < npoints; i++ )
  {
    f << std::setprecision(9) << x[i] << " " << y << " " << z << " "
      << velq[i](0) << " " << velq[i](1) << " " << velq[i](2) << " "
      << std::endl;
  }
  f.close();
  std::cout << "Wrote file source1.dat." << std::endl;

  f.open("source2.dat");
  f << "VARIABLES = \"x\" \"y\" \"z\" \"Vx\" \"Vy\" \"Vz\"" << std::endl;
  for ( i = 0; i < npoints; i++ )
  {
    f << std::setprecision(9) << x[i] << " " << y << " " << z << " "
      << velt[i](0) << " " << velt[i](1) << " " << velt[i](2) << " "
      << std::endl;
  }
  f.close();
  std::cout << "Wrote file source2.dat." << std::endl;

  f.open("sourcearr1.dat");
  f << "VARIABLES = \"x\" \"y\" \"z\" \"Vx\" \"Vy\" \"Vz\"" << std::endl;
  for ( i = 0; i < npoints; i++ )
  {
    f << std::setprecision(9) << x[i] << " " << y << " " << z << " "
      << velarr[i](0) << " " << velarr[i](1) << " " << velarr[i](2) << " "
      << std::endl;
  }
  f.close();
  std::cout << "Wrote file sourcearr1.dat." << std::endl;

  return 0;
}
