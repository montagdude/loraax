// Stores global settings

#define _USE_MATH_DEFINES

#include <string>
#include <Eigen/Core>
#include <cmath>
#include <tinyxml2.h>
extern "C"
{
  #include <xfoil_interface.h>
}
#include "util.h"
#include "settings.h"

using namespace tinyxml2;

// Case settings

std::string casename;
double uinf;
Eigen::Vector3d uinfvec;
double pinf;
double rhoinf;
double muinf;
double minf;
double alpha;
double dt;
double rollupdist;
int wakeiters;
double wakeangle;
bool viscous;
bool rollup_wake;
int reinit_freq;
double stop_tol;
int maxiters;
int viz_freq;

xfoil_geom_options_type xfoil_geom_opts;
xfoil_options_type xfoil_run_opts;

/******************************************************************************/
//
// Reads a single setting from XMLElement
//
/******************************************************************************/
int read_setting ( const XMLElement *elem, const std::string & setting,
                   std::string & value, bool required )
{
  if (elem->FirstChildElement(setting.c_str()))
    value = elem->FirstChildElement(setting.c_str())->GetText();
  else
  {
    if (required)
      conditional_stop(1, "read_setting",
                      "No " + setting +" in xml element " + elem->Name() + ".");
    return 1;
  }

  return 0;
}

int read_setting ( const XMLElement *elem, const std::string & setting,
                   double & value, bool required )
{
  if (elem->FirstChildElement(setting.c_str()))
    elem->FirstChildElement(setting.c_str())->QueryDoubleText(&value);
  else
  {
    if (required)
      conditional_stop(1, "read_setting",
                      "No " + setting +" in xml element " + elem->Name() + ".");
    return 1;
  }

  return 0;
}

int read_setting ( const XMLElement *elem, const std::string & setting,
                   int & value, bool required )
{
  if (elem->FirstChildElement(setting.c_str()))
    elem->FirstChildElement(setting.c_str())->QueryIntText(&value);
  else
  {
    if (required)
      conditional_stop(1, "read_setting",
                      "No " + setting +" in xml element " + elem->Name() + ".");
    return 1;
  }

  return 0;
}

int read_setting ( const XMLElement *elem, const std::string & setting,
                   bool & value, bool required )
{
  if (elem->FirstChildElement(setting.c_str()))
    elem->FirstChildElement(setting.c_str())->QueryBoolText(&value);
  else
  {
    if (required)
      conditional_stop(1, "read_setting",
                      "No " + setting +" in xml element " + elem->Name() + ".");
    return 1;
  }

  return 0;
}

/******************************************************************************/
//
// Read settings from input file
//
/******************************************************************************/
int read_settings ( const std::string & inputfile, std::string & geom_file )
{
	XMLDocument doc;

//FIXME: a lot of these settings need to be made optional with default values

	doc.LoadFile(inputfile.c_str());
	if ( (doc.ErrorID() == XML_ERROR_FILE_NOT_FOUND) ||
	     (doc.ErrorID() == XML_ERROR_FILE_COULD_NOT_BE_OPENED) ||
	     (doc.ErrorID() == XML_ERROR_FILE_READ_ERROR) )
	{
		conditional_stop(1, "read_settings",
		                 "Could not read " + inputfile + ".");
		return 1;
	}
	else if (doc.ErrorID() != 0)
	{
		conditional_stop(1, "read_settings",
		                 "Syntax error in " + inputfile + ".");
		return 1;
	}
	
	XMLElement *main = doc.FirstChildElement("Main");
	if (! main)
	{
		conditional_stop(1, "read_settings",
		                 "Expected 'Main' element in input file."); 
		return 2;
	}
	if (read_setting(main, "CaseName", casename) != 0)
		return 2;
	if (read_setting(main, "GeometryFile", geom_file) != 0)
		return 2;
	if (read_setting(main, "FreestreamSpeed", uinf) != 0)
		return 2;
	if (read_setting(main, "FreestreamStaticPressure", pinf) != 0)
		return 2;
	if (read_setting(main, "FreestreamDensity", rhoinf) != 0)
		return 2;
	if (read_setting(main, "FreestreamViscosity", muinf) != 0)
		return 2;
	if (read_setting(main, "AngleOfAttack", alpha) != 0)
		return 2;
	if (read_setting(main, "RollupWake", rollup_wake, false) != 0)
		rollup_wake = false;
	if (rollup_wake)
	{
		if (read_setting(main, "RollupDist", rollupdist) != 0)
			return 2;
		if (read_setting(main, "WakeIters", wakeiters) != 0)
			return 2;
	}
	else
	{
		rollupdist = -1.;	// Will get set to max span later
		wakeiters = 1;
	}
	if (read_setting(main, "InitialWakeAngle", wakeangle, false) != 0)
		wakeangle = alpha;
	if (read_setting(main, "Viscous", viscous) != 0)
		return 2;
	if (read_setting(main, "StoppingTolerance", stop_tol) != 0)
		return 2;
	if (read_setting(main, "MaxIters", maxiters) != 0)
		return 2;
	if (read_setting(main, "VisualizationFrequency", viz_freq) != 0)
		return 2;
	
	XMLElement *xfrun = main->FirstChildElement("XfoilRunOptions");
	if (! xfrun)
	{
		conditional_stop(1, "read_settings",
		                 "Expected 'XfoilRunOptions' element in input file.");
		return 2;
	}
	if (read_setting(xfrun, "ncrit", xfoil_run_opts.ncrit) != 0)
		return 2;
	if (read_setting(xfrun, "xtript", xfoil_run_opts.xtript) != 0)
		return 2;
	if (read_setting(xfrun, "xtripb", xfoil_run_opts.xtripb) != 0)
		return 2;
	if (read_setting(xfrun, "maxit", xfoil_run_opts.maxit) != 0)
		return 2;
	if (read_setting(xfrun, "vaccel", xfoil_run_opts.vaccel) != 0)
		return 2;
	if (read_setting(xfrun, "reinit_freq", reinit_freq, false) != 0)
		reinit_freq = 5;
	xfoil_run_opts.viscous_mode = viscous;
	xfoil_run_opts.silent_mode = true;
	
	XMLElement *xfgeom = main->FirstChildElement("XfoilPaneling");
	if (! xfgeom)
	{
		conditional_stop(1, "read_settings",
		                 "Expected 'XfoilPaneling' element in input file.");
		return 2;
	}
	if (read_setting(xfgeom, "npan", xfoil_geom_opts.npan) != 0)
		return 2;
	if (read_setting(xfgeom, "cvpar", xfoil_geom_opts.cvpar) != 0)
		return 2;
	if (read_setting(xfgeom, "cterat", xfoil_geom_opts.cterat) != 0)
		return 2;
	if (read_setting(xfgeom, "ctrrat", xfoil_geom_opts.ctrrat) != 0)
		return 2;
	xfoil_geom_opts.xsref1 = 1.0;
	xfoil_geom_opts.xsref2 = 1.0;
	xfoil_geom_opts.xpref1 = 1.0;
	xfoil_geom_opts.xpref2 = 1.0;
	
	// Set freestream vector, mach number, and time step size
	
	uinfvec(0) = uinf*cos(alpha*M_PI/180.);
	uinfvec(1) = 0.;
	uinfvec(2) = uinf*sin(alpha*M_PI/180.);
	minf = uinf / std::sqrt(1.4*pinf/rhoinf);
	if (rollup_wake)
		dt = rollupdist / (uinf * double(wakeiters));
	
	return 0;
}
