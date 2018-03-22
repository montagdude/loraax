#include <vector>
#include <Eigen/Core>
#include <tinyxml2.h>
#include "util.h"
#include "settings.h"
#include "section.h"
#include "wing.h"
#include "aircraft.h"

#include <iostream>

using namespace tinyxml2;

/******************************************************************************/
//
// Aircraft class. Contains some number of wings and related data and members.
// 
/******************************************************************************/

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
Aircraft::Aircraft ()
{
  _nwings = 0;
  _wings.resize(0);
  _sref = 0;
  _lref = 0;
  _momcen << 0., 0., 0.;
}

/******************************************************************************/
//
// Read from XML
//
/******************************************************************************/
int Aircraft::readXML ( const std::string & geom_file )
{
  XMLDocument doc;
  int nchord, nspan, check;
  double lesprat, tesprat, rootsprat, tipsprat;
  std::vector<Section> user_sections;
  std::vector<Airfoil> foils;
  double xle, y, zle, chord, twist, ymax;
  std::string source, des, path;
  const int npointside = 100;

  doc.LoadFile(geom_file.c_str());
  if ( (doc.ErrorID() == XML_ERROR_FILE_NOT_FOUND) ||
       (doc.ErrorID() == XML_ERROR_FILE_COULD_NOT_BE_OPENED) ||
       (doc.ErrorID() == XML_ERROR_FILE_READ_ERROR) )
  {
    conditional_stop(1, "Aircraft::readXML",
                     "Could not read " + geom_file + ".");
    return 1;
  }
  else if (doc.ErrorID() != 0)
  {
    conditional_stop(1, "Aircraft::readXML",
                     "Syntax error in " + geom_file + ".");
    return 1;
  }

  XMLElement *ac = doc.FirstChildElement("Aircraft");
  if (! ac)
  {
    conditional_stop(1, "Aircraft::readXML",
                     "Expected 'Aircraft' element in input file.");
    return 2;
  }

  // Read reference data

  XMLElement *ref = ac->FirstChildElement("Reference");
  if (! ac)
  {
    conditional_stop(1, "Aircraft::readXML",
                     "Expected 'Reference' element in input file.");
    return 2;
  }
  if (read_setting(ref, "RefArea", _sref) != 0)
    return 2;
  if (read_setting(ref, "RefLength", _lref) != 0)
    return 2;
  if (read_setting(ref, "MomentCenX", _momcen[0]) != 0)
    return 2;
  if (read_setting(ref, "MomentCenY", _momcen[1]) != 0)
    return 2;
  if (read_setting(ref, "MomentCenZ", _momcen[2]) != 0)
    return 2;

  // Read wing data

  for ( XMLElement *wingelem = ac->FirstChildElement("Wing"); wingelem != NULL;
        wingelem = wingelem->NextSiblingElement("Wing") )
  {
    Wing newwing;
    const char *wingname = wingelem->Attribute("name");
    if (wingname)
      newwing.setName(wingname);

    // Paneling

    XMLElement *pan = wingelem->FirstChildElement("Paneling");
    if (! pan)
    {
      conditional_stop(1, "Aircraft::readXML",
                       "Wing lacks 'Paneling' element.");
      return 2;
    }
    if (read_setting(pan, "NChord", nchord) != 0)
      return 2;
    if (read_setting(pan, "NSpan", nspan) != 0)
      return 2;
    if (read_setting(pan, "LESpaceRat", lesprat) != 0)
      return 2;
    if (read_setting(pan, "TESpaceRat", tesprat) != 0)
      return 2;
    if (read_setting(pan, "RootSpaceRat", rootsprat) != 0)
      return 2;
    if (read_setting(pan, "TipSpaceRat", tipsprat) != 0)
      return 2;
    newwing.setDiscretization(nchord, nspan, lesprat, tesprat, rootsprat,
                              tipsprat);

    // Sections

    user_sections.resize(0);
    ymax = 0.;
    XMLElement *secs = wingelem->FirstChildElement("Sections");
    if (! wingelem)
    {
      conditional_stop(1, "Aircraft::readXML",
                       "Wing lacks 'Sections' element.");
      return 2;
    }
    for ( XMLElement *secelem = secs->FirstChildElement("Section");
          secelem != NULL; secelem = secelem->NextSiblingElement("Section") )
    {
      Section newsection;
      const char *secname = secelem->Attribute("name");
      if ( (secname) && (std::string(secname) == "Root") )
        y = 0.;
      else
      {
        if (read_setting(secelem, "Y", y) != 0)
          return 2;
      }
      if (y < 0.)
      {
        conditional_stop(1, "Aircraft::readXML",
            "Y must be >= 0 for all sections. Mirroring occurs automatically.");
        return 2;
      }
      else if (y > ymax)
        ymax = y;
      if (read_setting(secelem, "XLE", xle) != 0)
        return 2;
      if (read_setting(secelem, "ZLE", zle) != 0)
        return 2;
      if (read_setting(secelem, "Chord", chord) != 0)
        return 2;
      if (read_setting(secelem, "Twist", twist) != 0)
        return 2;

      newsection.setGeometry(xle, y, zle, chord, twist, 0.0);
      user_sections.push_back(newsection);
    }
    if (user_sections.size() < 2)
    {
      conditional_stop(1, "Aircraft::readXML",
                       "Wings must have at least two sections.");
      return 2;
    }

    // Airfoils

    foils.resize(0);
    XMLElement *foilselem = wingelem->FirstChildElement("Airfoils");
    if (! foilselem)
    {
      conditional_stop(1, "Aircraft::readXML",
                       "Wing lacks 'Airfoils' element.");
      return 2;
    }
    for ( XMLElement *foilelem = foilselem->FirstChildElement("Airfoil");
          foilelem != NULL; foilelem = foilelem->NextSiblingElement("Airfoil") )
    {
      Airfoil newfoil;
      const char *foilname = foilelem->Attribute("name");
      if ( (foilname) && (std::string(foilname) == "Root") )
        y = 0.;
      else if ( (foilname) && (std::string(foilname) == "Tip") )
        y = ymax;
      else
      {
        if (read_setting(foilelem, "Y", y) != 0)
          return 2;
      }
      newfoil.setY(y);
      if (read_setting(foilelem, "Source", source) != 0)
        return 2;
      if (source == "4 digit")
      {
        if (read_setting(foilelem, "Designation", des) != 0)
          return 2;
        if (newfoil.naca4Coordinates(des, npointside) != 0)
        {
          conditional_stop(1, "Aircraft::readXML",
                           "Invalid 4-digit NACA designation.");
          return 2;
        }
      }
      else if (source == "5 digit")
      {
        if (read_setting(foilelem, "Designation", des) != 0)
          return 2;
        if (newfoil.naca5Coordinates(des, npointside) != 0)
        {
          conditional_stop(1, "Aircraft::readXML",
                           "Invalid 5-digit NACA designation.");
          return 2;
        }
      }
      else if (source == "file")
      {
        if (read_setting(foilelem, "Path", path) != 0)
          return 2;
        check = newfoil.readCoordinates(path);
        if (check == 1)
        {
          conditional_stop(1, "Aircraft::readXML",
                           "Error reading file " + path + ".");
          return 2;
        }
        else if (check == 2)
        {
          conditional_stop(1, "Aircraft::readXML",
                           "Format error in " + path + ".");
          return 2;
        }
      }
      // Set up airfoil spline data and smooth paneling

      newfoil.ccwOrderCoordinates();
      newfoil.splineFit();
      newfoil.unitTransform();
      newfoil.smoothPaneling(xfoil_geom_opts);
      foils.push_back(newfoil); 
    }
    newwing.setAirfoils(foils);    
    newwing.setupSections(user_sections);
    addWing(newwing); 
  }

  if (_nwings < 1)
  {
    conditional_stop(1, "Aircraft::readXML", "At least one wing is required.");
    return 2;
  }

  return 0;
}

/******************************************************************************/
//
// Adds a new wing
//
/******************************************************************************/
void Aircraft::addWing ( const Wing & wing )
{
  _nwings += 1;
  _wings.push_back(wing);
}
