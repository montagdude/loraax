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
#include "farfield.h"
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
    _allwake.resize(nwings);
    
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
        
        _allwake[i] = &_wings[i].wake();
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
  f.setf(std::ios_base::scientific);
  for ( i = 0; i < int(nverts); i++ )
  {
    f << std::setprecision(7) << std::setw(16) << std::left
      << _verts[i]->xViz();
    f << std::setprecision(7) << std::setw(16) << std::left
      << _verts[i]->yViz();
    f << std::setprecision(7) << std::setw(16) << std::left
      << _verts[i]->zViz() << std::endl;
  }
  for ( i = 0; i < int(nverts); i++ )
  {
    f << std::setprecision(7) << std::setw(16) << std::left
      << _verts[i]->xViz();
    f << std::setprecision(7) << std::setw(16) << std::left
      << -_verts[i]->yViz();
    f << std::setprecision(7) << std::setw(16) << std::left
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
    writeSurfaceScalar(f, "source_strength", 0);
    writeSurfaceScalar(f, "doublet_strength", 1);
    
    f << "Vectors velocity double" << std::endl;
    f.setf(std::ios_base::scientific);
    for ( i = 0; i < nverts; i++ )
    {
        f << std::setprecision(7) << std::setw(16) << std::left
          << _verts[i]->data(2);
        f << std::setprecision(7) << std::setw(16) << std::left
          << _verts[i]->data(3);
        f << std::setprecision(7) << std::setw(16) << std::left
          << _verts[i]->data(4) << std::endl;
    }
    for ( i = 0; i < nverts; i++ )
    {
        f << std::setprecision(7) << std::setw(16) << std::left
          << _verts[i]->data(2);
        f << std::setprecision(7) << std::setw(16) << std::left
          << -_verts[i]->data(3);
        f << std::setprecision(7) << std::setw(16) << std::left
          << _verts[i]->data(4) << std::endl;
    }
    
    writeSurfaceScalar(f, "pressure", 5);
    writeSurfaceScalar(f, "pressure_coefficient", 6);
    writeSurfaceScalar(f, "mach", 7);
    writeSurfaceScalar(f, "density", 8);
    if (viscous)
    {
        writeSurfaceScalar(f, "skin_friction_coefficient", 9);
        writeSurfaceScalar(f, "displacement_thickness", 10);
        writeSurfaceScalar(f, "log_amplification_ratio", 11);
        writeSurfaceScalar(f, "uedge", 12);
        writeSurfaceScalar(f, "cp2d", 13);
    }

    int j, npanels;
    npanels = _panels.size();
    f << "CELL_DATA " << npanels*2 << std::endl;
    f << "Scalars mass_defect_derivative double 1" << std::endl;
    f << "LOOKUP_TABLE default" << std::endl;
    for ( j = npanels-1; j >= 0; j-- )
    {
        f << std::setprecision(7) << _panels[j]->massDefectDerivative()
          << std::endl;
    }
    for ( j = 0; j < npanels; j++ )
    {
        f << std::setprecision(7) << _panels[j]->massDefectDerivative()
          << std::endl;
    }
}

/******************************************************************************/
//
// Writes surface scalar to VTK viz file
//
/******************************************************************************/
void Aircraft::writeSurfaceScalar ( std::ofstream & f,
                                    const std::string & varname,
                                    unsigned int varidx ) const
{
  unsigned int i, nverts;

  nverts = _verts.size();

  f << "SCALARS " << varname << " double 1" << std::endl;
  f << "LOOKUP_TABLE default" << std::endl;
  f.setf(std::ios_base::scientific);
  for ( i = 0; i < nverts; i++ )
  {
    f << std::setprecision(7) << _verts[i]->data(varidx) << std::endl;
  }
  for ( i = 0; i < nverts; i++ )
  {
    f << std::setprecision(7) << _verts[i]->data(varidx) << std::endl;
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
  f.setf(std::ios_base::scientific);
  for ( i = 0; i < int(nverts); i++ )
  {
    f << std::setprecision(7) << std::setw(16) << std::left
      << _wakeverts[i]->xViz();
    f << std::setprecision(7) << std::setw(16) << std::left
      << _wakeverts[i]->yViz();
    f << std::setprecision(7) << std::setw(16) << std::left
      << _wakeverts[i]->zViz() << std::endl;
  }
  for ( i = 0; i < int(nverts); i++ )
  {
    f << std::setprecision(7) << std::setw(16) << std::left
      << _wakeverts[i]->xViz();
    f << std::setprecision(7) << std::setw(16) << std::left
      << -_wakeverts[i]->yViz();
    f << std::setprecision(7) << std::setw(16) << std::left
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
    f.setf(std::ios_base::scientific);
    for ( i = 0; i < nwakeverts; i++ )
    {
        f << std::setprecision(7) << _wakeverts[i]->data(1) << std::endl;
    }
    for ( i = 0; i < nwakeverts; i++ )
    {
        f << std::setprecision(7) << _wakeverts[i]->data(1) << std::endl;
    }

    f << "SCALARS wake_time double 1" << std::endl;
    f << "LOOKUP_TABLE default" << std::endl;
    for ( i = 0; i < nwakeverts; i++ )
    {
        f << std::setprecision(7) << _wakeverts[i]->wakeTime() << std::endl;
    }
    for ( i = 0; i < nwakeverts; i++ )
    {
        f << std::setprecision(7) << _wakeverts[i]->wakeTime() << std::endl;
    }
}

/******************************************************************************/
//
// Writes legacy VTK farfield viz
//
/******************************************************************************/
int Aircraft::writeFarfieldViz ( const std::string & prefix )
{
  std::string ffname, fname;
  std::ofstream f;
  int i, j;
  unsigned int nverts, npanels, verts_offset, cellsize, ncellverts;

  ffname = prefix + "_farfield.vtk";
  fname = "visualization/" + ffname;

  f.open(fname.c_str());
  if (! f.is_open())
  {
    conditional_stop(1, "Aircraft::writeFarfieldViz",
                     "Unable to open " + fname + " for writing.");
    return 1;
  }

  // Header

  f << "# vtk DataFile Version 3.0" << std::endl;
  f << casename << std::endl;
  f << "ASCII" << std::endl;
  f << "DATASET UNSTRUCTURED_GRID" << std::endl;

  // Write vertices

  nverts = _farfield.nVerts();
  verts_offset = _verts.size() + _wakeverts.size();
  f << "POINTS " << nverts << " double" << std::endl;
  f.setf(std::ios_base::scientific);
  for ( i = 0; i < int(nverts); i++ )
  {
    f << std::setprecision(7) << std::setw(16) << std::left
      << _farfield.vert(i)->xViz();
    f << std::setprecision(7) << std::setw(16) << std::left
      << _farfield.vert(i)->yViz();
    f << std::setprecision(7) << std::setw(16) << std::left
      << _farfield.vert(i)->zViz() << std::endl;
  }

  // Write panels

  npanels = _farfield.nQuads();
  cellsize = 0;
  for ( i = 0; i < int(npanels); i++ )
  {
    cellsize += 1 + _farfield.quadPanel(i)->nVertices();
  }
  f << "CELLS " << npanels << " " << cellsize << std::endl;
  for ( i = 0; i < int(npanels); i++ )
  {
    ncellverts = _farfield.quadPanel(i)->nVertices();
    if ( (ncellverts != 4) && (ncellverts != 3) )
    {
      conditional_stop(1, "Aircraft::writeFarfieldViz",
                       "Farfield panels must all be quad- or tri-type.");
      return 1;
    }
    f << ncellverts;
    for ( j = 0; j < int(ncellverts); j++ )
    {
      f << " " << _farfield.quadPanel(i)->vertex(j).idx()-verts_offset;
    }
    f << std::endl;
  }

  // Cell types

  f << "CELL_TYPES " << npanels << std::endl;
  for ( i = 0; i < int(npanels); i++ )
  {
    ncellverts = _farfield.quadPanel(i)->nVertices();
    if (ncellverts == 4)
      f << 9 << std::endl;
    else if (ncellverts == 3)
      f << 5 << std::endl;
  }

  writeFarfieldData(f);
  f.close();

  return 0;
}

/******************************************************************************/
//
// Writes farfield data to VTK viz file
//
/******************************************************************************/
void Aircraft::writeFarfieldData ( std::ofstream & f )
{
    unsigned int i, nverts;
    
    nverts = _farfield.nVerts();
    
    // Vertex data
    
    f << "POINT_DATA " << nverts << std::endl;
    f << "Vectors velocity double" << std::endl;
    f.setf(std::ios_base::scientific);
    for ( i = 0; i < nverts; i++ )
    {
        f << std::setprecision(7) << std::setw(16) << std::left
          << _farfield.vert(i)->data(2);
        f << std::setprecision(7) << std::setw(16) << std::left
          << _farfield.vert(i)->data(3);
        f << std::setprecision(7) << std::setw(16) << std::left
          << _farfield.vert(i)->data(4) << std::endl;
    }
    writeFarfieldScalar(f, "pressure", 5);
    writeFarfieldScalar(f, "pressure_coefficient", 6);
    writeFarfieldScalar(f, "mach", 7);
    writeFarfieldScalar(f, "density", 8);
}

/******************************************************************************/
//
// Writes farfield scalar to VTK viz file
//
/******************************************************************************/
void Aircraft::writeFarfieldScalar ( std::ofstream & f,
                                     const std::string & varname,
                                     unsigned int varidx )
{
  unsigned int i, nverts;

  nverts = _farfield.nVerts();

  f << "SCALARS " << varname << " double 1" << std::endl;
  f << "LOOKUP_TABLE default" << std::endl;
  f.setf(std::ios_base::scientific);
  for ( i = 0; i < nverts; i++ )
  {
    f << std::setprecision(7) << _farfield.vert(i)->data(varidx) << std::endl;
  }
}

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
    _xte = 0.;
    _zte = 0.;
    _maxspan = 0.;
    _verts.resize(0);
    _panels.resize(0);
    _wakeverts.resize(0);
    _wakepanels.resize(0);
    _allwake.resize(0);
    _sourceic.resize(0,0);
    _doubletic.resize(0,0);
    _aic.resize(0,0);
    _mun.resize(0);
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
    double camber, xcamber, thick;
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
    for ( XMLElement *wingelem = ac->FirstChildElement("Wing");
          wingelem != NULL; wingelem = wingelem->NextSiblingElement("Wing") )
    {
        nwings += 1;
    }
    _wings.resize(nwings);
    
    // Read wing data
    
    nwings = 0;
    for ( XMLElement *wingelem = ac->FirstChildElement("Wing");
          wingelem != NULL; wingelem = wingelem->NextSiblingElement("Wing") )
    {
      const char *wingname = wingelem->Attribute("name");
        if (wingname)
            _wings[nwings].setName(wingname);
        else
            _wings[nwings].setName("Wing" + int2string(nwings+1));
        
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
              secelem != NULL;
              secelem = secelem->NextSiblingElement("Section") )
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
                    "Y must be >= 0 for all sections. Mirroring is automatic.");
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
              foilelem != NULL;
              foilelem = foilelem->NextSiblingElement("Airfoil") )
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
                    "Y must be >= 0 for all airfoils. Mirroring is automatic.");
                return 2;
            }
            newfoil.setY(y);
            if (read_setting(foilelem, "Source", source) != 0)
                return 2;
            if (source == "4 digit")
            {
                if (read_setting(foilelem, "Camber", camber) != 0)
                    return 2;
                if (read_setting(foilelem, "XCamber", xcamber) != 0)
                    return 2;
                if (read_setting(foilelem, "Thickness", thick) != 0)
                    return 2;
                if (newfoil.naca4Coordinates(camber, xcamber, thick,
                                             npointside) != 0)
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
            newfoil.setXfoilOptions(xfoil_run_opts, xfoil_geom_opts);
            newfoil.smoothPaneling();
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
    
    // Determine furthest aft root TE points for Trefftz plane calculation

    _xte = -1.E+12;
    _zte = 0.;
    _maxspan = 0.;
    for ( i = 0; i < nwings; i++ )
    {
        if (_wings[i].xTE() >= _xte)
            _xte = _wings[i].xTE();
        _zte += _wings[i].zTE();
        if (_wings[i].span() >= _maxspan)
            _maxspan = _wings[i].span();
    }
    _zte /= double(nwings);

    if (rollupdist < 0.)
    {
        rollupdist = _maxspan;
        dt = rollupdist / (uinf * double(wakeiters));
    }

    // Set up wake for each wing
    
    for ( i = 0; i < nwings; i++ )
    {
        _wings[i].setupWake(_maxspan, next_global_vertidx, next_global_elemidx,
                            i);;
    }
    
    if (nwings < 1)
    {
        conditional_stop(1, "Aircraft::readXML",
                         "At least one wing is required.");
        return 2;
    }

    // Farfield box

    if (enable_farfield)
        _farfield.initialize(farfield_nx, farfield_ny, farfield_nz,
                             farfield_cenx, farfield_ceny, farfield_cenz,
                             farfield_lenx, farfield_leny, farfield_lenz, minf,
                             next_global_vertidx, next_global_elemidx);
    
    // Set pointers to vertices, panels, and wake elements
    
    setGeometryPointers();
    
    return 0;
}

/******************************************************************************/
//
// Sets source strengths
//
/******************************************************************************/
void Aircraft::setSourceStrengths ( bool init )
{
    unsigned int i, npanels, nwings;

    npanels = _panels.size();
#ifdef DEBUG
    if (npanels == 0)
        conditional_stop(1, "Aircraft::setSourceStrengths", "No panels exist.");
#endif

#pragma omp parallel for private(i)
    for ( i = 0; i < npanels; i++ )
    {
        _panels[i]->computeSourceStrength(uinfvec, viscous);
    }

    // Update source strengths in viscous wake

    if ( (viscous) && ! (init) )
    {
        nwings = _wings.size();
        for ( i = 0; i < nwings; i++ )
        {
            _wings[i].viscousWake().update();
        }
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
    unsigned int i, j, k, npanels, nwings, nstrips, nwakepans, nwakeverts;
    double mu;
    WakeStrip *strip;

    npanels = _panels.size();
#ifdef DEBUG
    if (npanels == 0)
        conditional_stop(1, "Aircraft::setDoubletStrengths",
                         "No panels exist.");
    if (_mun.size() != npanels)
        conditional_stop(1, "Aircraft::setDoubletStrengths",
                     "Inconsistent number of panels and solution vector size.");
#endif

#pragma omp parallel for private(i)
    for ( i = 0; i < npanels; i++ )
    {
        _panels[i]->setDoubletStrength(_mun(i)*uinf);
    }

    // Set wake doublet strength

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
        _wakeverts[i]->averageFromPanels();
    }
}

/******************************************************************************/
//
// Constructs AIC matrix and RHS vector
//
/******************************************************************************/
void Aircraft::constructSystem ( bool init )
{
    unsigned int i, j, k, l, m, nwings, npanels, nstrips, nwakepans, nvwtris;
    int toptepan, bottepan;
    Eigen::Vector3d col;
    WakeStrip * strip;
    Panel * vwtri;
    double stripic;
    bool onpanel;

    npanels = _panels.size();
#ifdef DEBUG
    if (npanels == 0)
        conditional_stop(1, "Aircraft::constructSystem", "No panels exist.");
#endif
    nwings = _wings.size();

    // Compute influence coefficient matrices the first time through

    if (init)
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
                _sourceic(i,j) = _panels[j]->sourcePhiCoeff(col(0), col(1),
                                               col(2), onpanel, "bottom", true);
                _doubletic(i,j) = _panels[j]->doubletPhiCoeff(col(0), col(1),
                                               col(2), onpanel, "bottom", true);
            }
        }
    }

    // Compute AIC and RHS

#pragma omp parallel for private(i,col,j,k,nstrips,l,strip,nwakepans,stripic,\
                                 m,toptepan,bottepan,nvwtris,vwtri)
    for ( i = 0; i < npanels; i++ )
    {
        // Collocation point (point of BC application)
        
        col = _panels[i]->collocationPoint();
        
        // Surface panel contribution to AIC and RHS
        
        _rhs(i) = 0.;
        for ( j = 0; j < npanels; j++ )
        {
            _rhs(i) -= _panels[j]->sourceStrength()*_sourceic(i,j);
            if (init || rollup_wake)
                _aic(i,j) = _doubletic(i,j);
        }
        
        // Wake contribution to AIC and RHS
        
        for ( k = 0; k < nwings; k++ )
        {
            if (init || rollup_wake)
            {
                nstrips = _wings[k].nWStrips();
                for ( l = 0; l < nstrips; l++ )
                {
                    strip = _wings[k].wStrip(l);
                    nwakepans = strip->nPanels();
                    
                    // All wake panels in a strip have strength equal to
                    // mu_topte - mu_botte
                    
                    stripic = 0.;
                    for ( m = 0; m < nwakepans; m++ )
                    {
                        stripic += strip->panel(m)->doubletPhiCoeff(col(0),
                                         col(1), col(2), false, "bottom", true);
                    }
                    toptepan = strip->topTEPan()->idx();
                    bottepan = strip->botTEPan()->idx();
                    _aic(i,toptepan) += stripic;
                    _aic(i,bottepan) -= stripic;
                }
            }

            // Viscous wake influence

            if ( (viscous) && (! init) )
            {
                nvwtris = _wings[k].viscousWake().nTris();
                for ( l = 0; l < nvwtris; l++ )
                {
                    vwtri = _wings[k].viscousWake().triPanel(l);
                    _rhs(i) -= vwtri->sourcePhiCoeff(col(0), col(1), col(2),
                                  false, "bottom", true)
                             * vwtri->sourceStrength();
                }
            }
        }

        // Normalize RHS by uinf to keep magnitudes small in the linear system.
        // The resulting doublet strengths are later scaled back up.
        
        _rhs(i) /= uinf;
    }
}

/******************************************************************************/
//
// Factorizes the AIC matrix and solves the system
//
/******************************************************************************/
void Aircraft::factorize () { _lu.compute(_aic); }
void Aircraft::solveSystem () { _mun = _lu.solve(_rhs); }

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
    Eigen::Vector3d wi;
    
    nwings = _wings.size();
    for ( i = 0; i < nwings; i++ )
    {
        _wings[i].computeSurfaceQuantities();
    }
}

/******************************************************************************/
//
// BL calculations with xfoil
//
/******************************************************************************/
void Aircraft::computeBL ()
{
  unsigned int i, nwings;

  nwings = _wings.size();
  for ( i = 0; i < nwings; i++ )
  {
    _wings[i].computeBL();
  }
}

/******************************************************************************/
//
// Sets up viscous wake for each wing
//
/******************************************************************************/
void Aircraft::setupViscousWake ()
{
    unsigned int i, nwings;
    int next_global_vertidx, next_global_elemidx;

    next_global_vertidx = _verts.size() + _wakeverts.size();
    next_global_elemidx = _panels.size() + _wakepanels.size();
    nwings = _wings.size();
    for ( i = 0; i < nwings; i++ )
    {
        _wings[i].setupViscousWake(next_global_vertidx, next_global_elemidx);
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
    double beta, xteinc, zteinc, xtrefftz, ztrefftz;

    // P-G transformation to get aft TE point in incompressible coordinates

    beta = std::sqrt(1. - std::pow(minf,2.));
    xteinc = _xte/beta;
    zteinc = _zte;

    // Trefftz plane 1000 span lengths downstream from aftmost TE

    xtrefftz = xteinc + 1000.*_maxspan*uinfvec(0)/uinf;
    ztrefftz = zteinc + 1000.*_maxspan*uinfvec(2)/uinf;
    
    nwings = _wings.size();
    for ( i = 0; i < nwings; i++ )
    {
        _wings[i].computeForceMoment(_momcen, xtrefftz, ztrefftz, _allwake);
    }
}

double Aircraft::lift () const
{
    unsigned int i, nwings;
    double lift;

    nwings = _wings.size();
    lift = 0.;
    for ( i = 0; i < nwings; i++ )
    {
        lift += _wings[i].lift();
    }

    return lift;
}
double Aircraft::trefftzLift () const
{
    unsigned int i, nwings;
    double lifttr;

    nwings = _wings.size();
    lifttr = 0.;
    for ( i = 0; i < nwings; i++ )
    {
        lifttr += _wings[i].trefftzLift();
    }

    return lifttr;
}
double Aircraft::skinFrictionLift () const
{
    unsigned int i, nwings;
    double liftf;

    nwings = _wings.size();
    liftf = 0.;
    for ( i = 0; i < nwings; i++ )
    {
        liftf += _wings[i].skinFrictionLift();
    }

    return liftf;
}
double Aircraft::integratedLift () const
{
    unsigned int i, nwings;
    double liftp;

    nwings = _wings.size();
    liftp = 0.;
    for ( i = 0; i < nwings; i++ )
    {
        liftp += _wings[i].integratedLift();
    }

    return liftp;
}

double Aircraft::drag () const
{
    unsigned int i, nwings;
    double drag;

    nwings = _wings.size();
    drag = 0.;
    for ( i = 0; i < nwings; i++ )
    {
        drag += _wings[i].drag();
    }

    return drag;
}
double Aircraft::inducedDrag () const
{
    unsigned int i, nwings;
    double dragtr;

    nwings = _wings.size();
    dragtr = 0.;
    for ( i = 0; i < nwings; i++ )
    {
        dragtr += _wings[i].inducedDrag();
    }

    return dragtr;
}
double Aircraft::parasiticDrag () const
{
    unsigned int i, nwings;
    double dragv;

    nwings = _wings.size();
    dragv = 0.;
    for ( i = 0; i < nwings; i++ )
    {
        dragv += _wings[i].parasiticDrag();
    }

    return dragv;
}
double Aircraft::skinFrictionDrag () const
{
    unsigned int i, nwings;
    double dragf;

    nwings = _wings.size();
    dragf = 0.;
    for ( i = 0; i < nwings; i++ )
    {
        dragf += _wings[i].skinFrictionDrag();
    }

    return dragf;
}
double Aircraft::integratedDrag () const
{
    unsigned int i, nwings;
    double dragp;

    nwings = _wings.size();
    dragp = 0.;
    for ( i = 0; i < nwings; i++ )
    {
        dragp += _wings[i].integratedDrag();
    }

    return dragp;
}

double Aircraft::pitchingMoment () const
{
    unsigned int i, nwings;
    double moment;

    nwings = _wings.size();
    moment = 0.;
    for ( i = 0; i < nwings; i++ )
    {
        moment += _wings[i].pitchingMoment();
    }

    return moment;
}
double Aircraft::pressurePitchingMoment () const
{
    unsigned int i, nwings;
    double momentp;

    nwings = _wings.size();
    momentp = 0.;
    for ( i = 0; i < nwings; i++ )
    {
        momentp += _wings[i].pressurePitchingMoment();
    }

    return momentp;
}
double Aircraft::skinFrictionPitchingMoment () const
{
    unsigned int i, nwings;
    double momentf;

    nwings = _wings.size();
    momentf = 0.;
    for ( i = 0; i < nwings; i++ )
    {
        momentf += _wings[i].skinFrictionPitchingMoment();
    }

    return momentf;
}

double Aircraft::liftCoefficient () const
{
    double qinf;

    qinf = 0.5*rhoinf*std::pow(uinf, 2.);
    return lift() / (qinf*_sref);
}
double Aircraft::trefftzLiftCoefficient () const
{
    double qinf;

    qinf = 0.5*rhoinf*std::pow(uinf, 2.);
    return trefftzLift() / (qinf*_sref);
}
double Aircraft::skinFrictionLiftCoefficient () const
{
    double qinf;

    qinf = 0.5*rhoinf*std::pow(uinf, 2.);
    return skinFrictionLift() / (qinf*_sref);
}
double Aircraft::integratedLiftCoefficient () const
{
    double qinf;

    qinf = 0.5*rhoinf*std::pow(uinf, 2.);
    return integratedLift() / (qinf*_sref);
}

double Aircraft::dragCoefficient () const
{
    double qinf;

    qinf = 0.5*rhoinf*std::pow(uinf, 2.);
    return drag() / (qinf*_sref);
}
double Aircraft::inducedDragCoefficient () const
{
    double qinf;

    qinf = 0.5*rhoinf*std::pow(uinf, 2.);
    return inducedDrag() / (qinf*_sref);
}
double Aircraft::parasiticDragCoefficient () const
{
    double qinf;

    qinf = 0.5*rhoinf*std::pow(uinf, 2.);
    return parasiticDrag() / (qinf*_sref);
}
double Aircraft::skinFrictionDragCoefficient () const
{
    double qinf;

    qinf = 0.5*rhoinf*std::pow(uinf, 2.);
    return skinFrictionDrag() / (qinf*_sref);
}
double Aircraft::integratedDragCoefficient () const
{
    double qinf;

    qinf = 0.5*rhoinf*std::pow(uinf, 2.);
    return integratedDrag() / (qinf*_sref);
}

double Aircraft::pitchingMomentCoefficient () const
{
    double qinf;

    qinf = 0.5*rhoinf*std::pow(uinf, 2.);
    return pitchingMoment() / (qinf*_sref*_lref);
}
double Aircraft::pressurePitchingMomentCoefficient () const
{
    double qinf;

    qinf = 0.5*rhoinf*std::pow(uinf, 2.);
    return pressurePitchingMoment() / (qinf*_sref*_lref);
}
double Aircraft::skinFrictionPitchingMomentCoefficient () const
{
    double qinf;

    qinf = 0.5*rhoinf*std::pow(uinf, 2.);
    return skinFrictionPitchingMoment() / (qinf*_sref*_lref);
}

/******************************************************************************/
//
// Writes forces and moments for aircraft and each wing to CSV formatted files
//
/******************************************************************************/
int Aircraft::writeForceMoment ( int iter ) const
{
    std::ofstream f;
    std::string fname;
    unsigned int i, nwings;
    
    fname = "forcemoment/" + casename + "_forcemoment.csv";
    
    // Write header during first iteration
    
    if (iter == 1)
    {
        f.open(fname.c_str(), std::fstream::out);
        if (! f.is_open())
        {
            print_warning("Aircraft::writeForceMoment",
                          "Unable to open " + fname + " for writing.");
            return 1;
        }
        if (viscous)
        {
        f << "\"Iter\",\"Lift\",\"Lift_trefftz\",\"Lift_skinfric\","
          <<            "\"Drag\",\"Drag_induced\",\"Drag_parasitic\","
          <<"\"Moment\",\"Moment_pressure\",\"Moment_skinfric\","
          <<"\"CL\",\"CL_trefftz\",\"CL_skinfric\","
          <<"\"CD\",\"CD_induced\",\"CD_parasitic\","
          <<"\"Cm\",\"Cm_pressure\",\"Cm_skinfric\"" << std::endl;
        }
        else
        {
            f << "\"Iter\",\"Lift\",\"Lift_integrated\","
              <<          "\"Drag\",\"Drag_integrated\","
              << "\"Moment\","
              << "\"CL\",\"CL_integrated\","
              << "\"CD\",\"CD_integrated\","
              << "\"Cm\"" << std::endl;
        }
    }
    else
    {
        f.open(fname.c_str(), std::fstream::app);
        if (! f.is_open())
        {
          print_warning("Aircraft::writeForceMoment",
                        "Unable to open " + fname + " for writing.");
          return 1;
        }
    }
    f << iter << ",";
    f.setf(std::ios_base::scientific);
    f << std::setprecision(7);
    if (viscous)
    {
        f << lift() << ",";
        f << trefftzLift() << ",";
        f << skinFrictionLift() << ",";
        f << drag() << ",";
        f << inducedDrag() << ",";
        f << parasiticDrag() << ",";
        f << pitchingMoment() << ",";
        f << pressurePitchingMoment() << ",";
        f << skinFrictionPitchingMoment() << ",";
        f << liftCoefficient() << ",";
        f << trefftzLiftCoefficient() << ",";
        f << skinFrictionLiftCoefficient() << ",";
        f << dragCoefficient() << ",";
        f << inducedDragCoefficient() << ",";
        f << parasiticDragCoefficient() << ",";
        f << pitchingMomentCoefficient() << ",";
        f << pressurePitchingMomentCoefficient() << ",";
        f << skinFrictionPitchingMomentCoefficient() << std::endl;
    }
    else
    {
        f << lift() << ",";
        f << integratedLift() << ",";
        f << drag() << ",";
        f << integratedDrag() << ",";
        f << pitchingMoment() << ",";
        f << liftCoefficient() << ",";
        f << integratedLiftCoefficient() << ",";
        f << dragCoefficient() << ",";
        f << integratedDragCoefficient() << ",";
        f << pitchingMomentCoefficient() << std::endl;
    }
    f.close();
    
    nwings = _wings.size();
    for ( i = 0; i < nwings; i++ )
    {
        if (_wings[i].writeForceMoment(iter) != 0)
            return 1;
    }
    
    return 0;
}

/******************************************************************************/
//
// Writes sectional force and moment coefficients to file
//
/******************************************************************************/
void Aircraft::writeSectionForceMoment ( int iter ) const
{
  unsigned int i, nwings;

  nwings = _wings.size();
  for ( i = 0; i < nwings; i++ )
  {
    _wings[i].writeSectionForceMoment(iter);
  }
}

/******************************************************************************/
//
// Convects and updates wake panels
//
/******************************************************************************/
void Aircraft::moveWake ()
{
    unsigned int i, nwings;
    
    nwings = _wings.size();
    for ( i = 0; i < nwings; i++ )
    {
        _wings[i].wake().convectVertices(dt, _panels, _wakepanels);
    }
    for ( i = 0; i < nwings; i++ )
    {
        _wings[i].wake().update();
    }
}

/*******************************************************************************

Computes farfield data

*******************************************************************************/
void Aircraft::computeFarfield ()
{
  _farfield.computeVelocity(uinfvec, minf, _panels, _wakepanels);
  _farfield.computePressure(uinf, rhoinf, pinf);
  _farfield.computeForce(alpha, rhoinf, uinf, _sref);
}

/******************************************************************************/
//
// Writes legacy VTK viz files
//
/******************************************************************************/
int Aircraft::writeViz ( const std::string & prefix, int iter ) const
{
  std::string surfname, wakename;

  surfname = prefix + "_surfs_iter" + int2string(iter) + ".vtk";
  wakename = prefix + "_wake_iter" + int2string(iter) + ".vtk";

  if (writeSurfaceViz("visualization/" + surfname) != 0)
    return 1;

  if (writeWakeViz("visualization/" + wakename) != 0)
    return 1;

  return 0;
}

/*******************************************************************************

Writes farfield data

*******************************************************************************/
int Aircraft::writeFarfieldData ( const std::string & prefix ) const
{
    return _farfield.writeForceAccel(prefix);
}
