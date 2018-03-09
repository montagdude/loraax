#include <iostream>
extern "C"
{
  #include <xfoil_interface.h>
}
#include "settings.h"
#include "airfoil.h"

int main (int argc, char* argv[]) 
{ 
  Airfoil foil1, foil2, foil3;
  int npointside, check;
  xfoil_geom_options_type geom_opts;

  ncrit = 9.;
  xtript = 1.;
  xtripb = 1.;
  maxit = 100;
  vaccel = 0.01;
  fix_unconverged = true;
  reinitialize = true;

  npan = 160;
  cvpar = 1.;
  cterat = 0.15;
  ctrrat = 0.2;

  geom_opts.npan = npan;
  geom_opts.cvpar = cvpar;
  geom_opts.cterat = cterat;
  geom_opts.ctrrat = ctrrat;
  geom_opts.xsref1 = 1.;
  geom_opts.xsref2 = 1.;
  geom_opts.xpref1 = 1.;
  geom_opts.xpref2 = 1.;

  farfield_distance_factor = 1000.;

  check = foil1.readCoordinates("clarky.dat");
  if (check == 1)
    std::cout << "Unable to open clarky.dat." << std::endl;
  else if (check == 2)
    std::cout << "Format error in clarky.dat." << std::endl;
  foil1.ccwOrderCoordinates();
  foil1.unitTransform();
  foil1.smoothPaneling(geom_opts);

  npointside = 100;
  check = foil2.naca5Coordinates("21021", npointside);
  foil2.ccwOrderCoordinates();
  foil2.unitTransform();
  foil2.smoothPaneling(geom_opts);

  foil3.interpCoordinates(foil1, foil2, 0.5);
  foil3.unitTransform();
  foil3.smoothPaneling(geom_opts);
  std::cout << "Number of buffer points for foil3: " << foil3.nBuffer()
            << std::endl;

  return 0; 
}
