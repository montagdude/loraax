#define USE_MATH_DEFINES

#include <iomanip>
#include <vector>
#include <fstream>
#include <cmath>
#include <tinyxml2.h>
#include <Eigen/Core>
#include "util.h"
#include "settings.h"
#include "section.h"
#include "wake_strip.h"
#include "wake.h"
#include "element.h"
#include "panel.h"
#include "wing.h"
#include "aircraft.h"

using namespace tinyxml2;

/******************************************************************************/
//
// Aircraft class. Contains some number of wings and related data and members.
// 
/******************************************************************************/

/******************************************************************************/
//
// Sets up pointers to vertices, panels, and wake elements
//
/******************************************************************************/
void Aircraft::setGeometryPointers ()
{
  unsigned int nwings, nverts_total, nquads_total, ntris_total;
  unsigned int i, j, nverts, nquads, ntris, vcounter, pcounter;
  unsigned int nverts_wake_total, nverts_wake, vwcounter;
  unsigned int nwaketris_total, nwakequads_total, nwaketris, nwakequads;
  unsigned int wakecounter;

  // Get sizes first (don't use push_back, because it invalidates pointers)

  nwings = _wings.size();
  nverts_total = 0;
  nquads_total = 0;
  ntris_total = 0;
  nverts_wake_total = 0;
  nwaketris_total = 0;
  nwakequads_total = 0;
  for ( i = 0; i < nwings; i++ )
  {
    nverts_total += _wings[i].nVerts();
    nquads_total += _wings[i].nQuads();
    ntris_total += _wings[i].nTris();
    nverts_wake_total += _wings[i].wake().nVerts();
    nwaketris_total += _wings[i].wake().nTris();
    nwakequads_total += _wings[i].wake().nQuads();
  }
  _verts.resize(nverts_total);
  _panels.resize(nquads_total + ntris_total);
  _wakeverts.resize(nverts_wake_total);
  _wakepanels.resize(nwaketris_total + nwakequads_total);

  // Store geometry pointers

  vcounter = 0;
  pcounter = 0;
  vwcounter = 0;
  wakecounter = 0;
  for ( i = 0; i < nwings; i++ ) 
  {
    nverts = _wings[i].nVerts();
    for ( j = 0; j < nverts; j++ )
    {
      _verts[vcounter] = _wings[i].vert(j);
      vcounter += 1;
    }

    nquads = _wings[i].nQuads();
    for ( j = 0; j < nquads; j++ )
    {
      _panels[pcounter] = _wings[i].quadPanel(j);
      pcounter += 1;
    }

    ntris = _wings[i].nTris();
    for ( j = 0; j < ntris; j++ )
    {
      _panels[pcounter] = _wings[i].triPanel(j);
      pcounter += 1;
    }

    nverts_wake = _wings[i].wake().nVerts();
    for ( j = 0; j < nverts_wake; j++ )
    {
      _wakeverts[vwcounter] = _wings[i].wake().vert(j);
      vwcounter += 1;
    }

    nwaketris = _wings[i].wake().nTris();
    for ( j = 0; j < nwaketris; j++ )
    {
      _wakepanels[wakecounter] = _wings[i].wake().triPanel(j);
      wakecounter += 1;
    }

    nwakequads = _wings[i].wake().nQuads();
    for ( j = 0; j < nwakequads; j++ )
    {
      _wakepanels[wakecounter] = _wings[i].wake().quadPanel(j);
      wakecounter += 1;
    } 
  }
}

/******************************************************************************/
//
// Writes legacy VTK surface viz
//
/******************************************************************************/
int Aircraft::writeSurfaceViz ( const std::string & fname ) const
{
  std::ofstream f;
  int i, j;
  unsigned int nverts, npanels, cellsize, ncellverts;

  f.open(fname.c_str());
  if (! f.is_open())
  {
    conditional_stop(1, "Aircraft::writeSurfaceViz",
                     "Unable to open " + fname + " for writing."); 
    return 1;
  }

  // Header

  f << "# vtk DataFile Version 3.0" << std::endl;
  f << casename << std::endl;
  f << "ASCII" << std::endl;
  f << "DATASET UNSTRUCTURED_GRID" << std::endl;

  // Write vertices and mirror vertices

  nverts = _verts.size();
  f << "POINTS " << nverts*2 << " double" << std::endl;
  for ( i = 0; i < int(nverts); i++ )
  {
    f << std::setprecision(14) << std::setw(25) << std::left
      << _verts[i]->xViz();
    f << std::setprecision(14) << std::setw(25) << std::left
      << _verts[i]->yViz();
    f << std::setprecision(14) << std::setw(25) << std::left
      << _verts[i]->zViz() << std::endl;
  } 
  for ( i = 0; i < int(nverts); i++ )
  {
    f << std::setprecision(14) << std::setw(25) << std::left
      << _verts[i]->xViz();
    f << std::setprecision(14) << std::setw(25) << std::left
      << -_verts[i]->yViz();
    f << std::setprecision(14) << std::setw(25) << std::left
      << _verts[i]->zViz() << std::endl;
  } 

  // Write panels and mirror panels

  npanels = _panels.size();
  cellsize = 0;
  for ( i = 0; i < int(npanels); i++ )
  {
    cellsize += 1 + _panels[i]->nVertices();
  }
  f << "CELLS " << npanels*2 << " " << cellsize*2 << std::endl;
  for ( i = npanels-1; i >= 0; i-- )
  {
    ncellverts = _panels[i]->nVertices();    
    if ( (ncellverts != 4) && (ncellverts != 3) )
    {
      conditional_stop(1, "Aircraft::writeSurfaceViz",
                       "Panels must all be quad- or tri-type.");  
      return 1;
    }
    f << ncellverts;
    for ( j = 0; j < int(ncellverts); j++ )
    {
      f << " " << _panels[i]->vertex(j).idx();
    }
    f << std::endl;
  }
  for ( i = 0; i < int(npanels); i++ )
  {
    ncellverts = _panels[i]->nVertices();    
    f << ncellverts;
    for ( j = ncellverts-1; j >= 0; j-- )
    {
      f << " " << _panels[i]->vertex(j).idx()+nverts;
    }
    f << std::endl;
  }

  // Cell types for panels and mirror panels

  f << "CELL_TYPES " << npanels*2 << std::endl;
  for ( i = npanels-1; i >= 0; i-- )
  {
    ncellverts = _panels[i]->nVertices();    
    if (ncellverts == 4)
      f << 9 << std::endl;
    else if (ncellverts == 3)
      f << 5 << std::endl;
  }
  for ( i = 0; i < int(npanels); i++ )
  {
    ncellverts = _panels[i]->nVertices();    
    if (ncellverts == 4)
      f << 9 << std::endl;
    else if (ncellverts == 3)
      f << 5 << std::endl;
  }

  // Surface data at vertices

  writeSurfaceData(f);

  f.close();

  return 0;
}

/******************************************************************************/
//
// Writes surface data to VTK viz file
//
/******************************************************************************/
void Aircraft::writeSurfaceData ( std::ofstream & f ) const
{
  unsigned int i, nverts;

  nverts = _verts.size();
  
  // Vertex data (incl. mirror panels)

  f << "POINT_DATA " << nverts*2 << std::endl;
  f << "SCALARS source_strength double 1" << std::endl;
  f << "LOOKUP_TABLE default" << std::endl;
  for ( i = 0; i < nverts; i++ )
  {
    f << std::setprecision(14) << std::setw(25) << std::left
      << _verts[i]->data(0) << std::endl;
  }
  for ( i = 0; i < nverts; i++ )
  {
    f << std::setprecision(14) << std::setw(25) << std::left
      << _verts[i]->data(0) << std::endl;
  }

  f << "SCALARS doublet_strength double 1" << std::endl;
  f << "LOOKUP_TABLE default" << std::endl;
  for ( i = 0; i < nverts; i++ )
  {
    f << std::setprecision(14) << std::setw(25) << std::left
      << _verts[i]->data(1) << std::endl;
  }
  for ( i = 0; i < nverts; i++ )
  {
    f << std::setprecision(14) << std::setw(25) << std::left
      << _verts[i]->data(1) << std::endl;
  }

  f << "Vectors velocity double" << std::endl;
  for ( i = 0; i < nverts; i++ )
  {
    f << std::setprecision(14) << std::setw(25) << std::left
      << _verts[i]->data(2);
    f << std::setprecision(14) << std::setw(25) << std::left
      << _verts[i]->data(3);
    f << std::setprecision(14) << std::setw(25) << std::left
      << _verts[i]->data(4) << std::endl;
  }
  for ( i = 0; i < nverts; i++ )
  {
    f << std::setprecision(14) << std::setw(25) << std::left
      << _verts[i]->data(2);
    f << std::setprecision(14) << std::setw(25) << std::left
      << -_verts[i]->data(3);
    f << std::setprecision(14) << std::setw(25) << std::left
      << _verts[i]->data(4) << std::endl;
  }

  f << "SCALARS pressure double 1" << std::endl;
  f << "LOOKUP_TABLE default" << std::endl;
  for ( i = 0; i < nverts; i++ )
  {
    f << std::setprecision(14) << std::setw(25) << std::left
      << _verts[i]->data(5) << std::endl;
  }
  for ( i = 0; i < nverts; i++ )
  {
    f << std::setprecision(14) << std::setw(25) << std::left
      << _verts[i]->data(5) << std::endl;
  }

  f << "SCALARS pressure_coefficient double 1" << std::endl;
  f << "LOOKUP_TABLE default" << std::endl;
  for ( i = 0; i < nverts; i++ )
  {
    f << std::setprecision(14) << std::setw(25) << std::left
      << _verts[i]->data(6) << std::endl;
  }
  for ( i = 0; i < nverts; i++ )
  {
    f << std::setprecision(14) << std::setw(25) << std::left
      << _verts[i]->data(6) << std::endl;
  }
}

/******************************************************************************/
//
// Writes legacy VTK wake viz
//
/******************************************************************************/
int Aircraft::writeWakeViz ( const std::string & fname ) const
{
  std::ofstream f;
  int i, j;
  unsigned int nverts, npanels, nsurf_verts, cellsize, ncellverts;

  f.open(fname.c_str());
  if (! f.is_open())
  {
    conditional_stop(1, "Aircraft::writeWakeViz",
                     "Unable to open " + fname + " for writing."); 
    return 1;
  }

  // Header

  f << "# vtk DataFile Version 3.0" << std::endl;
  f << casename << std::endl;
  f << "ASCII" << std::endl;
  f << "DATASET UNSTRUCTURED_GRID" << std::endl;

  // Write vertices and mirror vertices

  nverts = _wakeverts.size();
  nsurf_verts = _verts.size();
  f << "POINTS " << nverts*2 << " double" << std::endl;
  for ( i = 0; i < int(nverts); i++ )
  {
    f << std::setprecision(14) << std::setw(25) << std::left
      << _wakeverts[i]->xViz();
    f << std::setprecision(14) << std::setw(25) << std::left
      << _wakeverts[i]->yViz();
    f << std::setprecision(14) << std::setw(25) << std::left
      << _wakeverts[i]->zViz() << std::endl;
  } 
  for ( i = 0; i < int(nverts); i++ )
  {
    f << std::setprecision(14) << std::setw(25) << std::left
      << _wakeverts[i]->xViz();
    f << std::setprecision(14) << std::setw(25) << std::left
      << -_wakeverts[i]->yViz();
    f << std::setprecision(14) << std::setw(25) << std::left
      << _wakeverts[i]->zViz() << std::endl;
  } 

  // Write wake panels and mirror panels

  npanels = _wakepanels.size();
  cellsize = 0;
  for ( i = 0; i < int(npanels); i++ )
  {
    cellsize += 1 + _wakepanels[i]->nVertices();
  }
  f << "CELLS " << npanels*2 << " " << cellsize*2 << std::endl;
  for ( i = npanels-1; i >= 0; i-- )
  {
    ncellverts = _wakepanels[i]->nVertices();
    if ( (ncellverts != 4) && (ncellverts != 3) )
    {
      conditional_stop(1, "Aircraft::writeWakeViz",
                       "Wake panels must all be quad- or tri-type.");
      return 1;
    }
    f << ncellverts;
    for ( j = 0; j < int(ncellverts); j++ )
    {
      f << " " << _wakepanels[i]->vertex(j).idx()-nsurf_verts;
    }
    f << std::endl;
  }
  for ( i = 0; i < int(npanels); i++ )
  {
    ncellverts = _wakepanels[i]->nVertices();
    f << ncellverts;
    for ( j = ncellverts-1; j >= 0; j-- )
    {
      f << " " << _wakepanels[i]->vertex(j).idx()-nsurf_verts+nverts;
    }
    f << std::endl;
  }

  // Cell types for doublet panels and mirror panels

  f << "CELL_TYPES " << npanels*2 << std::endl;
  for ( i = npanels-1; i >= 0; i-- )
  {
    ncellverts = _wakepanels[i]->nVertices();
    if (ncellverts == 4)
      f << 9 << std::endl;
    else if (ncellverts == 3) 
      f << 5 << std::endl;
  }
  for ( i = 0; i < int(npanels); i++ )
  {
    if (ncellverts == 4)
      f << 9 << std::endl;
    else if (ncellverts == 3) 
      f << 5 << std::endl;
  }

  // Wake data

  writeWakeData(f);

  f.close();

  return 0;
}

/******************************************************************************/
//
// Writes wake data to VTK viz file
//
/******************************************************************************/
void Aircraft::writeWakeData ( std::ofstream & f ) const
{
  unsigned int i, nwakeverts;

  nwakeverts = _wakeverts.size();
  
  // Point data (incl. mirror elements)

  f << "POINT_DATA " << nwakeverts*2 << std::endl;
  f << "SCALARS doublet_strength double 1" << std::endl;
  f << "LOOKUP_TABLE default" << std::endl;
  for ( i = 0; i < nwakeverts; i++ )
  {
    f << std::setprecision(14) << std::setw(25) << std::left
      << _wakeverts[i]->data(1) << std::endl;
  }
  for ( i = 0; i < nwakeverts; i++ )
  {
    f << std::setprecision(14) << std::setw(25) << std::left
      << _wakeverts[i]->data(1) << std::endl;
  }
}

#if 0
/******************************************************************************/
//
// Writes legacy VTK wake strip viz (for checking purposes). One file is written
// for each strip so they can be cycled through.
//
/******************************************************************************/
int Aircraft::writeWakeStripViz ( const std::string & prefix )
{
  std::ofstream f;
  unsigned int counter;
  unsigned int i, j, k, nwings, nstrips, nsurf_verts, nwake_verts, nvorts;
  std::string fname;
  std::vector<WakeStrip *> wakestrips;

  // Determine number of wake strips and populate container

  nwings = _wings.size();
  nstrips = 0;
  for ( i = 0; i < nwings; i++ )
  {
    nstrips += _wings[i].nWStrips();
  }
  wakestrips.resize(nstrips);
  counter = 0;
  for ( i = 0; i < nwings; i++ )
  {
    for ( j = 0; j < _wings[i].nWStrips(); j++ )
    {
      wakestrips[counter] = _wings[i].wStrip(j);
      counter += 1;
    }
  }

  // Write a file for each strip

  for ( j = 0; j < nstrips; j++ )
  {
    fname = prefix + "_timestep" + int2string(j) + ".vtk";
    f.open(fname.c_str());
    if (! f.is_open())
    {
      conditional_stop(1, "Aircraft::writeWakeStripViz",
                       "Unable to open " + fname + " for writing."); 
      return 1;
    }

    // Header

    f << "# vtk DataFile Version 3.0" << std::endl;
    f << casename << std::endl;
    f << "ASCII" << std::endl;
    f << "DATASET UNSTRUCTURED_GRID" << std::endl;

    // Write all surface and wake vertices

    nsurf_verts = _verts.size();
    nwake_verts = _wakeverts.size();
    f << "POINTS " << nsurf_verts + nwake_verts << " double" << std::endl;
    for ( i = 0; i < nsurf_verts; i++ )
    {
      f << std::setprecision(14) << std::setw(25) << std::left
        << _verts[i]->x();
      f << std::setprecision(14) << std::setw(25) << std::left
        << _verts[i]->y();
      f << std::setprecision(14) << std::setw(25) << std::left
        << _verts[i]->z() << std::endl;
    } 
    for ( i = 0; i < nwake_verts; i++ )
    {
      f << std::setprecision(14) << std::setw(25) << std::left
        << _wakeverts[i]->x();
      f << std::setprecision(14) << std::setw(25) << std::left
        << _wakeverts[i]->y();
      f << std::setprecision(14) << std::setw(25) << std::left
        << _wakeverts[i]->z() << std::endl;
    } 

    // Write TE panels and wake strip elements

    nvorts = wakestrips[j]->nVortices();
    f << "CELLS " << 2+nvorts << " " << 5*(2+nvorts) << std::endl;
    f << 4;
    for ( k = 0; k < 4; k++ )
    {
      f << " " << wakestrips[j]->topTEPan()->vertex(k).idx();
    }
    f << std::endl;
    f << 4;
    for ( k = 0; k < 4; k++ )
    {
      f << " " << wakestrips[j]->botTEPan()->vertex(k).idx();
    }
    f << std::endl;
    for ( i = 0; i < nvorts; i++ )
    {
      f << 4;
      for ( k = 0; k < 4; k++ )
      {
        f << " " << wakestrips[j]->vortex(i)->vertex(k).idx();
      }
      f << std::endl;
    } 

    // Cell types

    f << "CELL_TYPES " << 2+nvorts << std::endl;
    f << 9 << std::endl;
    f << 9 << std::endl;
    for ( i = 0; i < nvorts-1; i++ )
    {
      f << 9 << std::endl;
    }
    f << 4 << std::endl;

    f.close();
  }

  return 0;
}
#endif

/******************************************************************************/
//
// Default constructor
//
/******************************************************************************/
Aircraft::Aircraft ()
{
  _wings.resize(0);
  _sref = 0;
  _lref = 0;
  _momcen << 0., 0., 0.;
  _verts.resize(0);
  _panels.resize(0);
  _wakeverts.resize(0);
  _wakepanels.resize(0);
  _sourceic.resize(0,0);
  _doubletic.resize(0,0);
  _aic.resize(0,0);
  _mu.resize(0);
  _rhs.resize(0);
}

/******************************************************************************/
//
// Read from XML
//
/******************************************************************************/
int Aircraft::readXML ( const std::string & geom_file )
{
  XMLDocument doc;
  unsigned int i, nwings;
  int nchord, nspan, check;
  double lesprat, tesprat, rootsprat, tipsprat;
  std::vector<Section> user_sections;
  std::vector<Airfoil> foils;
  double xle, y, zle, chord, twist, ymax;
  std::string source, des, path;
  const int npointside = 100;
  int next_global_elemidx = 0;
  int next_global_vertidx = 0;

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

  // First determine number of wings and allocate _wings vector (don't use
  // push_back, because it invalidates pointers)

  nwings = 0;
  for ( XMLElement *wingelem = ac->FirstChildElement("Wing"); wingelem != NULL;
        wingelem = wingelem->NextSiblingElement("Wing") )
  {
    nwings += 1;
  }
  _wings.resize(nwings);

  // Read wing data

  nwings = 0;
  for ( XMLElement *wingelem = ac->FirstChildElement("Wing"); wingelem != NULL;
        wingelem = wingelem->NextSiblingElement("Wing") )
  {
    const char *wingname = wingelem->Attribute("name");
    if (wingname)
      _wings[nwings].setName(wingname);

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
    _wings[nwings].setDiscretization(nchord, nspan, lesprat, tesprat,
                                     rootsprat, tipsprat);

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

      // Some geometry checks

      if (chord <= 0.)
        conditional_stop(1, "Aircraft::readXML",
                         "Chord must be greater than 0.");
      if ( (twist <= -90.) || (twist >= 90.) )
        conditional_stop(1, "Aircraft::readXML",
                    "Twist must be between greater than -90 and less than 90.");

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
      if (y < 0)
      {
        conditional_stop(1, "Aircraft::readXML",
            "Y must be >= 0 for all airfoils. Mirroring occurs automatically.");
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
      else
      {
        conditional_stop(1, "Aircraft::readXML",
                         "Airfoil source must be 4 digit, 5 digit, or file.");
        return 2;
      }

      // Set up airfoil spline data and smooth paneling

      newfoil.ccwOrderCoordinates();
      newfoil.splineFit();
      newfoil.unitTransform();
      newfoil.smoothPaneling(xfoil_geom_opts);
      foils.push_back(newfoil); 
    }
    if (foils.size() < 1)
    {
      conditional_stop(1, "Aircraft::readXML",
                       "Wings must have at least one airfoil.");
       return 2;
    }

    _wings[nwings].setAirfoils(foils);    
    _wings[nwings].setupSections(user_sections);
    _wings[nwings].createPanels(next_global_vertidx, next_global_elemidx);
    nwings += 1;
  }

  // Set up wake for each wing

  for ( i = 0; i < nwings; i++ )
  {
    _wings[i].setupWake(next_global_vertidx, next_global_elemidx);
  }

  if (nwings < 1)
  {
    conditional_stop(1, "Aircraft::readXML", "At least one wing is required.");
    return 2;
  }

  // Set pointers to vertices, panels, and wake elements

  setGeometryPointers();

  return 0;
}

/******************************************************************************/
//
// Sets source strengths
//
/******************************************************************************/
void Aircraft::setSourceStrengths ()
{
  unsigned int i, npanels;
  Eigen::Vector3d norm;

  npanels = _panels.size();
#ifdef DEBUG
  if (npanels == 0)
    conditional_stop(1, "Aircraft::setSourceStrengths", "No panels exist.");
#endif

#pragma omp parallel for private(i,norm)
  for ( i = 0; i < npanels; i++ )
  {
    norm = _panels[i]->normal();
    _panels[i]->setSourceStrength(-uinfvec.transpose() * norm);
  }
}

/******************************************************************************/
//
// Sets doublet strengths from solution vector. Also interpolates doublet
// strengths to vertices.
//
/******************************************************************************/
void Aircraft::setDoubletStrengths ()
{
  unsigned int i, npanels;

  npanels = _panels.size();
#ifdef DEBUG
  if (npanels == 0)
    conditional_stop(1, "Aircraft::setDoubletStrengths", "No panels exist.");
  if (_mu.size() != npanels)
    conditional_stop(1, "Aircraft::setDoubletStrengths",
                     "Inconsistent number of panels and solution vector size.");
#endif

#pragma omp parallel for private(i)
  for ( i = 0; i < npanels; i++ )
  {
    _panels[i]->setDoubletStrength(_mu(i));
  }
}

/******************************************************************************/
//
// Sets wake doublet strength
//
/******************************************************************************/
void Aircraft::setWakeDoubletStrengths ()
{
  unsigned int i, j, k, nwings, nstrips, nwakepans, nwakeverts;
  double mu; 
  WakeStrip *strip;

  nwings = _wings.size();
  for ( i = 0; i < nwings; i++ )
  {
    nstrips = _wings[i].nWStrips();
#pragma omp parallel for private(j,strip,mu,nwakepans,k)
    for ( j = 0; j < nstrips; j++ )
    { 
      strip = _wings[i].wStrip(j);
      mu = strip->topTEPan()->doubletStrength()
         - strip->botTEPan()->doubletStrength();
      nwakepans = strip->nPanels();
      for ( k = 0; k < nwakepans; k++ )
      {
        strip->panel(k)->setDoubletStrength(mu);
      }
    }
  }

  nwakeverts = _wakeverts.size();
#pragma omp parallel for private(i)
  for ( i = 0; i < nwakeverts; i++ )
  {
    _wakeverts[i]->interpFromPanels();
  }
}

/******************************************************************************/
//
// Constructs AIC matrix and RHS vector
//
/******************************************************************************/
void Aircraft::constructSystem ()
{
  unsigned int i, j, k, l, m, nwings, npanels, nstrips, nwakepans;
  int toptepan, bottepan;
  Eigen::Vector3d col;
  WakeStrip * strip;
  double stripic; 
  bool onpanel;

  npanels = _panels.size();
#ifdef DEBUG
  if (npanels == 0)
    conditional_stop(1, "Aircraft::constructSystem", "No panels exist.");
#endif
  nwings = _wings.size();

  // Compute influence coefficient matrices the first time through

  if (_sourceic.rows() != npanels)
  {
    _sourceic.resize(npanels,npanels);
    _doubletic.resize(npanels,npanels);
    _aic.resize(npanels,npanels);
    _rhs.resize(npanels);

#pragma omp parallel for private(i,col,j,onpanel)
    for ( i = 0; i < npanels; i++ )
    {
      // Collocation point (point of BC application)

      col = _panels[i]->collocationPoint();

      // Influence coefficients

      for ( j = 0; j < npanels; j++ )
      {
        if (i == j)
          onpanel = true;
        else
          onpanel = false;
        _sourceic(i,j) = _panels[j]->sourcePhiCoeff(col(0), col(1), col(2),
                                                    onpanel, "bottom", true);
        _doubletic(i,j) = _panels[j]->doubletPhiCoeff(col(0), col(1), col(2),
                                                      onpanel, "bottom", true);
      }
    }
  }

  // Compute AIC and RHS

#pragma omp parallel for private(i,col,j,k,nstrips,l,strip,nwakepans,stripic,m,\
                                 toptepan,bottepan)
  for ( i = 0; i < npanels; i++ )
  {
    // Collocation point (point of BC application)

    col = _panels[i]->collocationPoint();

    // Surface panel contribution to AIC and RHS

    _rhs(i) = 0.;
    for ( j = 0; j < npanels; j++ )
    {
      _aic(i,j) = _doubletic(i,j);
      _rhs(i) -= _panels[j]->sourceStrength()*_sourceic(i,j);
    }

    // Wake contribution to AIC

    for ( k = 0; k < nwings; k++ )
    {
      nstrips = _wings[k].nWStrips(); 
      for ( l = 0; l < nstrips; l++ )
      {
        strip = _wings[k].wStrip(l);
        nwakepans = strip->nPanels();
        stripic = 0.;
        for ( m = 0; m < nwakepans; m++ )
        {
          stripic += strip->panel(m)->doubletPhiCoeff(col(0), col(1), col(2),
                                                      false, "bottom", true);
        }
        toptepan = strip->topTEPan()->idx();
        bottepan = strip->botTEPan()->idx();
        _aic(i,toptepan) += stripic;
        _aic(i,bottepan) -= stripic;
      }
    }
  }
}

/******************************************************************************/
//
// Factorizes the AIC matrix and solves the system
//
/******************************************************************************/
void Aircraft::factorize () { _lu.compute(_aic); }
void Aircraft::solveSystem () { _mu = _lu.solve(_rhs); }

/******************************************************************************/
//
// Gives size of system of equations (= number of panels)
//
/******************************************************************************/
unsigned int Aircraft::systemSize () const { return _panels.size(); }

/******************************************************************************/
//
// Computes surface velocities and pressures
//
/******************************************************************************/
void Aircraft::computeSurfaceQuantities ()
{
  unsigned int i, nwings;

  nwings = _wings.size();
  for ( i = 0; i < nwings; i++ )
  {
    _wings[i].computeSurfaceQuantities();
  }
}

/******************************************************************************/
//
// Computes or access forces and moments
//
/******************************************************************************/
void Aircraft::computeForceMoment ()
{
  unsigned int i, nwings;
  double qinf;

  nwings = _wings.size();
  _force << 0., 0., 0.;
  _moment << 0., 0., 0.;
  for ( i = 0; i < nwings; i++ )
  {
    _wings[i].computeForceMoment(_momcen);
    _force += _wings[i].force();
    _moment += _wings[i].moment();
  }

  // Wind frame forces and coefficients

  qinf = 0.5*rhoinf*std::pow(uinf, 2.);
  _lift = -_force(0)*sin(alpha*M_PI/180.) + _force(2)*cos(alpha*M_PI/180.);
  _drag =  _force(0)*cos(alpha*M_PI/180.) + _force(2)*sin(alpha*M_PI/180.);
  _cl = _lift/(qinf*_sref);
  _cd = _drag/(qinf*_sref);
  _cm = _moment(1)/(qinf*_sref*_lref);
}

const double & Aircraft::liftCoefficient () const { return _cl; }
const double & Aircraft::dragCoefficient () const { return _cd; }
const double & Aircraft::pitchingMomentCoefficient () const { return _cm; }

/******************************************************************************/
//
// Writes legacy VTK viz files
//
/******************************************************************************/
int Aircraft::writeViz ( const std::string & prefix ) const
{
  std::string surfname, wakename;

  surfname = prefix + "_surfs.vtk";
  wakename = prefix + "_wake.vtk";

  if (writeSurfaceViz(surfname) != 0)
    return 1;

  if (writeWakeViz(wakename) != 0)
    return 1;

  return 0;
}
