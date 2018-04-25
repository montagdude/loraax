#include <iostream>
#include <iomanip>
#include <string>
#include "clo_parser.h"
#include "settings.h"
#include "aircraft.h"

int main (int argc, char* argv[]) 
{ 
  const std::string version = "0.1";
  CLOParser parser;
  std::string geom_file;
  int check;
  Aircraft ac;
  unsigned int iter;

  // Parse CLOs

  check = parser.checkCLOs(argc, argv, version);
  if (check < 0)
    return 0;
  else if (check > 0)
    return 1;
  else
  {
    // Read settings

    std::cout << "Reading settings ..." << std::endl;
    if (read_settings(parser.inputFile(), geom_file) != 0)
      return 2;
  }
    
  // Read geometry

  std::cout << "Reading and discretizing geometry ..." << std::endl;
  if (ac.readXML(geom_file) != 0)
    return 3;

  iter = 1;
  while (int(iter) <= maxsteps)
  {
    // Convect wake

    std::cout << "Iteration " << iter << std::endl;

    if (iter > 1)
    {
      std::cout << "  Convecting wake ..." << std::endl;
      ac.moveWake();
    }

    // Set source strengths

    std::cout << "  Setting source strengths ..." << std::endl;
    ac.setSourceStrengths();

    // Construct, factorize, and solve the system

    std::cout << "  Constructing the linear system ..." << std::endl;
    ac.constructSystem(iter==1);
    std::cout << "  Factorizing the AIC matrix ..." << std::endl;
    ac.factorize();
    std::cout << "  Solving the linear system with " << ac.systemSize()
              << " unknowns ..." << std::endl;
    ac.solveSystem();

    // Set doublet strengths on surface and wake

    std::cout << "  Setting doublet strengths ..." << std::endl;
    ac.setDoubletStrengths();
    ac.setWakeDoubletStrengths(iter==1);

    // Compute surface velocities and pressures

    std::cout << "  Computing surface velocity and pressure ..." << std::endl;
    ac.computeSurfaceQuantities();

    // Compute forces and moments

    std::cout << "  Computing forces and moments ..." << std::endl;
    ac.computeForceMoment();
    std::cout << "  CL: " << std::setprecision(5) << std::setw(8) << std::left
              << ac.liftCoefficient();
    std::cout << "  CD: " << std::setprecision(5) << std::setw(8) << std::left
              << ac.dragCoefficient();
    std::cout << "  Cm: " << std::setprecision(5) << std::setw(8) << std::left
              << ac.pitchingMomentCoefficient() << std::endl;

    // Write viz
//FIXME: do this only every viz_freq iterations
    ac.writeViz(casename, iter);

    iter++;
  }

  return 0; 
}
