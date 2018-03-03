// Reads/writes data from/to files, including settings, grid, and visualization

#ifndef FILE_IO_H
#define FILE_IO_H

#include <string>
#include "body.h"
#include "wake.h"

// Public functions
  
void read_settings ( const std::string & ); 
                                         // Reads settings from namelist file

void read_grid ( const std::string &, const std::string &, Body & );
                                         // Reads geometry from grid file

void write_solution_output ( const std::string &, const std::string &,
                             const std::vector<std::string> &, Body &, Wake & );
                                         // Writes solution output to files

#endif
