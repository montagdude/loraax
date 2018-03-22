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

  std::cout << "Reading geometry ..." << std::endl;
  if (ac.readXML(geom_file) != 0)
    return 3;

  return 0; 
}
