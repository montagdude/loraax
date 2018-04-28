// Stores global settings

#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>
#include <Eigen/Core>
#include <tinyxml2.h>
extern "C"
{
  #include <xfoil_interface.h>
}

using namespace tinyxml2;

// Case settings

extern std::string casename;
extern double uinf;
extern Eigen::Vector3d uinfvec;
extern double pinf;
extern double rhoinf;
extern double minf;
extern double muinf;
extern double alpha;
extern double dt;
extern double rollupdist;
extern double wakeangle;
extern bool viscous;
extern bool compressible;
extern double stop_tol;
extern int maxsteps;
extern int viz_freq;

// Xfoil settings

extern xfoil_geom_options_type xfoil_geom_opts;
extern xfoil_options_type xfoil_run_opts;

// Functions

int read_setting ( const XMLElement *elem, const std::string & setting,
                   std::string & value, bool required=true );
int read_setting ( const XMLElement *elem, const std::string & setting,
                   double & value, bool required=true );
int read_setting ( const XMLElement *elem, const std::string & setting,
                   int & value, bool required=true );
int read_setting ( const XMLElement *elem, const std::string & setting,
                   bool & value, bool required=true );
int read_settings ( const std::string & inputfile, std::string & geom_file );

#endif
