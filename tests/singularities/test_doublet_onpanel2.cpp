#include <vector>
#include <Eigen/Core>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "vertex.h"
#include "quadpanel.h"
#include "settings.h"

int main ()
{
  Vertex v1, v2, v3, v4;
  QuadPanel quadpanel;
  std::vector<double> phi, z;
  unsigned int i, npoints;
  double x, y, zmin, zmax, dz;
  std::ofstream f;

  // Evaluates potential in a vertical line passing through the doublet to make
  // sure the on-panel value is correct

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

  // First approach from below

  npoints = 50;
  zmin = -3.;
  zmax = 0.;
  dz = (zmax-zmin)/double(npoints-1);
  z.resize(npoints*2);
  phi.resize(npoints*2);
  x = 0.5;
  y = -0.25;
  for ( i = 0; i < npoints; i++ )
  {
    z[i] = zmin + double(i)*dz;
    if (std::abs(z[i]) < 1.E-12)
      phi[i] = quadpanel.inducedPotential(x, y, z[i], true, "bottom", false);
    else
      phi[i] = quadpanel.inducedPotential(x, y, z[i], false, "bottom", false);
  } 

  // Now approach from above

  zmin = 0.;
  zmax = 3.;
  dz = (zmax-zmin)/double(npoints-1);
  for ( i = 0; i < npoints; i++ )
  {
    z[i+npoints] = zmin + double(i)*dz;
    if (std::abs(z[i+npoints]) < 1.E-12)
      phi[i+npoints] = quadpanel.inducedPotential(x, y, z[i+npoints], true,
                                                  "top", false);
    else
      phi[i+npoints] = quadpanel.inducedPotential(x, y, z[i+npoints], false,
                                                  "top", false);
  } 

  f.open("quad_doublet2_below.dat");
  f << "VARIABLES = \"x\" \"y\" \"z\" \"Phi\"" << std::endl;
  for ( i = 0; i < npoints; i++ )
  {
    f << std::setprecision(9) << x << " " << y << " " << z[i] << " "
      << phi[i] << std::endl;
  }
  f.close();
  std::cout << "Wrote file quad_doublet2_below.dat." << std::endl;

  f.open("quad_doublet2_above.dat");
  f << "VARIABLES = \"x\" \"y\" \"z\" \"Phi\"" << std::endl;
  for ( i = 0; i < npoints; i++ )
  {
    f << std::setprecision(9) << x << " " << y << " " << z[i+npoints] << " "
      << phi[i+npoints] << std::endl;
  }
  f.close();
  std::cout << "Wrote file quad_doublet2_above.dat." << std::endl;

  return 0;
}
