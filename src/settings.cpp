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
int miniters;
int viz_freq;

xfoil_geom_options_type xfoil_geom_opts;
xfoil_options_type xfoil_run_opts;

bool enable_farfield;
double farfield_cenx, farfield_ceny, farfield_cenz;
double farfield_lenx, farfield_leny, farfield_lenz;
int farfield_nx, farfield_ny, farfield_nz;

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
                          "No " + setting +" in xml element " + elem->Name() +
                          ".");
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
                            "No " + setting +" in xml element " + elem->Name() +
                            ".");
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
                            "No " + setting +" in xml element " + elem->Name() +
                            ".");
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
                            "No " + setting +" in xml element " + elem->Name() +
                            ".");
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
    wakeangle = alpha;
    if (rollup_wake)
    {
        if (read_setting(main, "RollupDist", rollupdist) != 0)
            return 2;
        if (read_setting(main, "WakeIters", wakeiters) != 0)
            return 2;
        read_setting(main, "InitialWakeAngle", wakeangle, false);
    }
    else
    {
        rollupdist = -1.;   // Will get set to max span later
        wakeiters = 1;
    }
    if (read_setting(main, "Viscous", viscous) != 0)
        return 2;
    if (read_setting(main, "StoppingTolerance", stop_tol, false) != 0)
        stop_tol = 1.E-5;
    if (read_setting(main, "MaxIters", maxiters, false) != 0)
        maxiters = 100;
    if (read_setting(main, "MinIters", miniters, false) != 0)
        miniters = 0;
    if (read_setting(main, "VisualizationFrequency", viz_freq, false) != 0)
        viz_freq = 1;
    
    xfoil_run_opts.ncrit = 9.;
    xfoil_run_opts.xtript = 1.0;
    xfoil_run_opts.xtripb = 1.0;
    xfoil_run_opts.maxit = 100;
    xfoil_run_opts.vaccel = 0.01;
    reinit_freq = 5;
    XMLElement *xfrun = main->FirstChildElement("XfoilRunOptions");
    if (xfrun)
    {
        read_setting(xfrun, "ncrit", xfoil_run_opts.ncrit, false);
        read_setting(xfrun, "xtript", xfoil_run_opts.xtript, false);
        read_setting(xfrun, "xtripb", xfoil_run_opts.xtripb, false);
        read_setting(xfrun, "maxit", xfoil_run_opts.maxit, false);
        read_setting(xfrun, "vaccel", xfoil_run_opts.vaccel, false);
        read_setting(xfrun, "reinit_freq", reinit_freq, false);
    }
    xfoil_run_opts.viscous_mode = viscous;
    xfoil_run_opts.silent_mode = true;
    
    xfoil_geom_opts.npan = 160.;
    xfoil_geom_opts.cvpar = 1.0;
    xfoil_geom_opts.cterat = 0.15;
    XMLElement *xfgeom = main->FirstChildElement("XfoilPaneling");
    if (xfgeom)
    {
        read_setting(xfgeom, "npan", xfoil_geom_opts.npan, false);
        read_setting(xfgeom, "cvpar", xfoil_geom_opts.cvpar, false);
        read_setting(xfgeom, "cterat", xfoil_geom_opts.cterat, false);
    }
    xfoil_geom_opts.ctrrat = 0.20;
    xfoil_geom_opts.xsref1 = 1.0;
    xfoil_geom_opts.xsref2 = 1.0;
    xfoil_geom_opts.xpref1 = 1.0;
    xfoil_geom_opts.xpref2 = 1.0;

    // Postprocessing settings

    enable_farfield = false;
    XMLElement *post = main->FirstChildElement("Postprocessing");
    if (post)
    {
        XMLElement *farfield = post->FirstChildElement("Farfield");
        if (farfield)
        {
            if (read_setting(farfield, "Enable", enable_farfield) != 0)
                return 2;
            if (enable_farfield)
            {
              if (read_setting(farfield, "CenX", farfield_cenx) != 0)
                  return 2;
              if (read_setting(farfield, "CenY", farfield_ceny) != 0)
                  return 2;
              if (read_setting(farfield, "CenZ", farfield_cenz) != 0)
                  return 2;
              if (read_setting(farfield, "LenX", farfield_lenx) != 0)
                  return 2;
              if (read_setting(farfield, "LenY", farfield_leny) != 0)
                  return 2;
              if (read_setting(farfield, "LenZ", farfield_lenz) != 0)
                  return 2;
              if (read_setting(farfield, "NPointsX", farfield_nx) != 0)
                  return 2;
              if (read_setting(farfield, "NPointsY", farfield_ny) != 0)
                  return 2;
              if (read_setting(farfield, "NPointsZ", farfield_nz) != 0)
                  return 2;
            }
        }
    }
    
    // Set freestream vector, mach number, and time step size
    
    uinfvec(0) = uinf*cos(alpha*M_PI/180.);
    uinfvec(1) = 0.;
    uinfvec(2) = uinf*sin(alpha*M_PI/180.);
    minf = uinf / std::sqrt(1.4*pinf/rhoinf);
    if (minf >= 1.)
        conditional_stop(1, "read_settings",
                         "Mach number must be subsonic.");
    if (rollup_wake)
        dt = rollupdist / (uinf * double(wakeiters));
    
    return 0;
}
