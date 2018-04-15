#include <iostream>
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

  // Set source strengths

  std::cout << "Setting source strengths ..." << std::endl;
  ac.setSourceStrengths();

  // Construct, factorize, and solve the system

  std::cout << "Constructing the linear system ..." << std::endl;
  ac.constructSystem();
  std::cout << "Factorizing the AIC matrix ..." << std::endl;
  ac.factorize();
  std::cout << "Solving the linear system with " << ac.systemSize()
            << " unknowns ..." << std::endl;
  ac.solveSystem();

  // Set doublet strengths on surface and wake

  std::cout << "Setting doublet strengths ..." << std::endl;
  ac.setDoubletStrengths();
  ac.setWakeDoubletStrengths();

#if 0
  // Compute panel velocities

  std::cout << "Computing velocities ..." << std::endl;
  ac.computeVelocities();

  // Compute vertex data

  std::cout << "Computing vertex data ..." << std::endl;
  ac.computeVertexData();
#endif

  // Write viz

  ac.writeViz(casename);

  return 0; 
}
