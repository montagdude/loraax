#include <iostream>
#include "settings.h"
#include "airfoil.h"

int main (int argc, char* argv[]) 
{ 
  Airfoil foil;
  int check;

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

  farfield_distance_factor = 1000.;

  check = foil.readCoordinates("clarky.dat");
  if (check == 1)
    std::cout << "Unable to open clarky.dat." << std::endl;
  else if (check == 2)
    std::cout << "Format error in clarky.dat." << std::endl;
  foil.ccOrderCoordinates();

  return 0; 
}
