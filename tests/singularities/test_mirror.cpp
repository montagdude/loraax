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
  Vertex v1, v2, v3, v4, v5, v6;
  TriPanel tri1, tri2;
  std::vector<double> x;
  std::vector<Eigen::Vector3d> velmirror, veltwo;
  unsigned int npoints, i;
  double y, z, xmin, xmax, dx;
  std::ofstream f;

  v1.setCoordinates(0., 0.5, 0.);
  v2.setCoordinates(1., 0.8, 0.2);
  v3.setCoordinates(1., 1.5, 0.9);

  v4.setCoordinates(0., -0.5, 0.);
  v5.setCoordinates(1., -0.8, 0.2);
  v6.setCoordinates(1., -1.5, 0.9);

  tri1.addVertex(&v1);
  tri1.addVertex(&v2);
  tri1.addVertex(&v3);
  tri1.setSourceStrength(1.0);
  tri1.setDoubletStrength(1.0);

  tri2.addVertex(&v6);
  tri2.addVertex(&v5);
  tri2.addVertex(&v4);
  tri2.setSourceStrength(1.0);
  tri2.setDoubletStrength(1.0);
 
  npoints = 100;
  xmin = -5;
  xmax = 5;
  dx = (xmax-xmin)/double(npoints-1);
  x.resize(npoints);
  velmirror.resize(npoints);
  veltwo.resize(npoints);
  y = 0.0;
  z = -0.1;
  for ( i = 0; i < npoints; i++ )
  {
    x[i] = xmin + double(i)*dx;
    velmirror[i] = tri1.inducedVelocity(x[i], y, z, false, "top", true);
    veltwo[i] = tri1.inducedVelocity(x[i], y, z, false, "top", false)
              + tri2.inducedVelocity(x[i], y, z, false, "top", false);
  } 

  f.open("mirror1.dat");
  f << "VARIABLES = \"x\" \"y\" \"z\" \"Vx\" \"Vy\" \"Vz\"" << std::endl;
  for ( i = 0; i < npoints; i++ )
  {
    f << std::setprecision(9) << x[i] << " " << y << " " << z << " "
      << velmirror[i](0) << " " << velmirror[i](1) << " " << velmirror[i](2)
      << std::endl;
  }
  f.close();
  std::cout << "Wrote file mirror1.dat." << std::endl;

  f.open("mirror2.dat");
  f << "VARIABLES = \"x\" \"y\" \"z\" \"Vx\" \"Vy\" \"Vz\"" << std::endl;
  for ( i = 0; i < npoints; i++ )
  {
    f << std::setprecision(9) << x[i] << " " << y << " " << z << " "
      << veltwo[i](0) << " " << veltwo[i](1) << " " << veltwo[i](2)
      << std::endl;
  }
  f.close();
  std::cout << "Wrote file mirror2.dat." << std::endl;

  return 0;
}
