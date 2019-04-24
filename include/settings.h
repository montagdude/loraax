// Stores global settings

#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>
#include <Eigen/Core>
#include <tinyxml2.h>
extern "C"
{
  #include <libxfoil.h>
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
extern int wakeiters;
extern double wakeangle;
extern bool viscous;
extern bool rollup_wake;
extern int reinit_freq;
extern double stop_tol;
extern int maxiters;
extern int miniters;
extern int viz_freq;

// Xfoil settings

extern xfoil_geom_options_type xfoil_geom_opts;
extern xfoil_options_type xfoil_run_opts;

// Farfield settings

extern bool enable_farfield;
extern double farfield_cenx, farfield_ceny, farfield_cenz;
extern double farfield_lenx, farfield_leny, farfield_lenz;
extern int farfield_nx, farfield_ny, farfield_nz;

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
