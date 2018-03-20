#include <iostream>
#include <string>
#include "clo_parser.h"
#include "settings.h"

int main (int argc, char* argv[]) 
{ 
  const std::string version = "0.1";
  CLOParser parser;
  std::string geom_file;
  int check;

  check = parser.checkCLOs(argc, argv, version);
  if (check < 0)
    return 0;
  else if (check > 0)
    return 1;
  else
    read_settings(parser.inputFile(), geom_file);

std::cout << casename << std::endl;
std::cout << uinf << std::endl;
std::cout << rhoinf << std::endl;
std::cout << muinf << std::endl;
std::cout << alpha << std::endl;
std::cout << dt << std::endl;
std::cout << wakelen << std::endl;
std::cout << viscous << std::endl;
std::cout << stop_tol << std::endl;
std::cout << maxsteps << std::endl;
std::cout << viz_freq << std::endl;
std::cout << xfoil_run_opts.ncrit << std::endl;
std::cout << xfoil_run_opts.xtript << std::endl;
std::cout << xfoil_run_opts.xtripb << std::endl;
std::cout << xfoil_run_opts.maxit << std::endl;
std::cout << xfoil_run_opts.vaccel << std::endl;
std::cout << xfoil_geom_opts.npan << std::endl;
std::cout << xfoil_geom_opts.cvpar << std::endl;
std::cout << xfoil_geom_opts.cterat << std::endl;
std::cout << xfoil_geom_opts.ctrrat << std::endl;

  return 0; 
}
