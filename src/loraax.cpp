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

  return 0; 
}
