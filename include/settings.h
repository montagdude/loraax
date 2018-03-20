// Stores global settings

#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>
extern "C"
{
  #include <xfoil_interface.h>
}

// Case settings

extern std::string casename;
extern double uinf;
extern double rhoinf;
extern double muinf;
extern double alpha;
extern double dt;
extern double wakelen;
extern bool viscous;
extern double stop_tol;
extern int maxsteps;
extern int viz_freq;

// Xfoil settings

extern xfoil_geom_options_type xfoil_geom_opts;
extern xfoil_options_type xfoil_run_opts;

// Algorithm settings

extern double farfield_distance_factor;

// Functions

int read_settings ( const std::string & inputfile, std::string & geom_file );

#endif
