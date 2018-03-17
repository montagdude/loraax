#include <iostream>
extern "C"
{
  #include <xfoil_interface.h>
}
#include "settings.h"
#include "section.h"
#include "airfoil.h"

int main (int argc, char* argv[]) 
{ 
  Section sec1, sec2, sec3;
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

  check = sec1.airfoil().readCoordinates("clarky.dat");
  if (check == 1)
    std::cout << "Unable to open clarky.dat." << std::endl;
  else if (check == 2)
    std::cout << "Format error in clarky.dat." << std::endl;
  sec1.airfoil().ccwOrderCoordinates();
  sec1.airfoil().splineFit();
  sec1.airfoil().unitTransform();
  sec1.airfoil().smoothPaneling(geom_opts);

  npointside = 100;
  check = sec2.airfoil().naca5Coordinates("21021", npointside);
  sec2.airfoil().ccwOrderCoordinates();
  sec2.airfoil().splineFit();
  sec2.airfoil().unitTransform();
  sec2.airfoil().smoothPaneling(geom_opts);

  sec3.airfoil().interpCoordinates(sec1.airfoil(), sec2.airfoil(), 0.5);
  sec3.airfoil().splineFit();
  sec3.airfoil().unitTransform();
  sec3.airfoil().smoothPaneling(geom_opts);
  sec3.setVertices(51, 0.2, 0.5);

  return 0; 
}
