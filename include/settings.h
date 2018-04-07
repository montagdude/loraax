// Stores global settings

#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>
#include <tinyxml2.h>
extern "C"
{
  #include <xfoil_interface.h>
}

using namespace tinyxml2;

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

// Functions

int read_setting ( const XMLElement *elem, const std::string & setting,
                   std::string & value );
int read_setting ( const XMLElement *elem, const std::string & setting,
                   double & value );
int read_setting ( const XMLElement *elem, const std::string & setting,
                   int & value );
int read_setting ( const XMLElement *elem, const std::string & setting,
                   bool & value );
int read_settings ( const std::string & inputfile, std::string & geom_file );

#endif
