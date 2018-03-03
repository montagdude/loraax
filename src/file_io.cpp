// Reads/writes data from/to files, including settings, grid, and visualizations

#define _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <string>
#include <string.h>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include "namelist.h"
#include "settings.h"
#include "transformations.h"
#include "util.h"
#include "body.h"
#include "wake.h"

/******************************************************************************/
//
// Reads settings from namelist file
//
/******************************************************************************/
void read_settings ( const std::string & filename )
{
  std::string message, templist;
  std::ifstream infile;
  std::string wakelen;
  Namelist grid("grid", 3);
  Namelist oper("operating_conditions", 2);
  Namelist fomo("force_moment_reference", 7);
  Namelist grou("grouping", 2);
  Namelist wake("wake", 3);
  Namelist sout("solution_output", 4);
  Namelist algo("algorithms", 3);
  Eigen::Vector3d uinfw;
  Eigen::Matrix3d wind2body;

  // Open input file

  infile.open(filename.c_str());
  if (not infile.is_open())
  {
    message = "Input file " + filename + " could not be found.";
    conditional_stop(1, "read_settings", message);
  }

  // Set inputs and defaults (grid)

  grid.addItem("grid_format", "gmsh");
  grid.addItem("grid_file");
  grid.addItem("reverse_normals", "false");

  // Read namelist (grid)
  
  grid.readNamelist(infile);

  // Save settings (grid)

  grid_format = grid.getStringValue("grid_format");
  grid_file = grid.getStringValue("grid_file");
  reverse_normals = grid.getBoolValue("reverse_normals");

  // Set inputs and defaults (operating conditions)
    
  oper.addItem("alpha", "0.0");
  oper.addItem("beta", "0.0");

  // Read namelist (operating conditions)
  
  oper.readNamelist(infile);

  // Save settings (operating conditions)

  alpha = oper.getDoubleValue("alpha");
  beta = oper.getDoubleValue("beta");

  // Get freestream velocity in body frame from from alpha and beta

  uinfw(0) = 1.0;
  uinfw(1) = 0.0;
  uinfw(2) = 0.0;
  wind2body = euler_rotation(0.0, alpha*M_PI/180., beta*M_PI/180.);
  uinf = wind2body*uinfw;

  // Set inputs and defaults (forces and moments)

  fomo.addItem("ref_area", "1.0");
  fomo.addItem("mom_len_x", "1.0");
  fomo.addItem("mom_len_y", "1.0");
  fomo.addItem("mom_len_z", "1.0");
  fomo.addItem("mom_cen_x", "0.0");
  fomo.addItem("mom_cen_y", "0.0");
  fomo.addItem("mom_cen_z", "0.0");

  // Read namelist (forces and moments)

  fomo.readNamelist(infile);

  // Save settings (forces and moments)
  
  ref_area = fomo.getDoubleValue("ref_area");
  mom_len_x = fomo.getDoubleValue("mom_len_x");
  mom_len_y = fomo.getDoubleValue("mom_len_y");
  mom_len_z = fomo.getDoubleValue("mom_len_z");
  mom_cen_x = fomo.getDoubleValue("mom_cen_x");
  mom_cen_y = fomo.getDoubleValue("mom_cen_y");
  mom_cen_z = fomo.getDoubleValue("mom_cen_z");

  // Set inputs and defaults (grouping)

  grou.addItem("trailing_edge_list");
  grou.addItem("thin_surface_list");

  // Read namelist (grouping)
  
  grou.readNamelist(infile);

  // Save settings (grouping)

  templist = grou.getStringValue("trailing_edge_list");
  trailing_edge_list = split_string(templist, ',');
  templist = grou.getStringValue("thin_surface_list");
  thin_surface_list = split_string(templist, ',');

  // Set inputs and defaults (wake)

  wake.addItem("deforming_wake_length", double2string(10.*mom_len_y));
  wake.addItem("num_filaments", "10");
  wake.addItem("core_radius_factor", "0.25");

  // Read namelist (wake)
  
  wake.readNamelist(infile);

  // Save settings (wake)

  deforming_wake_length = wake.getDoubleValue("deforming_wake_length");
  num_filaments = wake.getIntValue("num_filaments");
  core_radius_factor = wake.getDoubleValue("core_radius_factor");

  // Set inputs and defaults (solution output)

  sout.addItem("case_name", "default_case");
  sout.addItem("output_file_format", "ensight");
  sout.addItem("output_variable_list");
  sout.addItem("write_edges", "false");

  // Read namelist (solution output)

  sout.readNamelist(infile);

  // Save settings (solution output)
  
  case_name = sout.getStringValue("case_name");
  output_file_format = sout.getStringValue("output_file_format");
  templist = sout.getStringValue("output_variable_list");
  output_variable_list = split_string(templist, ',');
  write_edges = sout.getBoolValue("write_edges");

  // Set inputs and defaults (algorithms)

  algo.addItem("farfield_distance_factor", "5.0");
  algo.addItem("relax_wake", "true");
  algo.addItem("wake_convergence_tolerance", "1.E-06");

  // Read namelist (algorithm)

  algo.readNamelist(infile);

  // Save settings (algorithms)

  farfield_distance_factor = algo.getDoubleValue("farfield_distance_factor");
  relax_wake = algo.getBoolValue("relax_wake");
  wake_convergence_tolerance = algo.getDoubleValue(
                                                  "wake_convergence_tolerance");

  // Close input file
  
  infile.close();

  // Echo namelist settings

  std::cout << "Echoing program settings:" << std::endl;
  std::cout << std::endl;
  grid.echo();
  std::cout << std::endl;
  oper.echo();
  std::cout << std::endl;
  fomo.echo();
  std::cout << std::endl;
  grou.echo();
  std::cout << std::endl;
  wake.echo();
  std::cout << std::endl;
  sout.echo();
  std::cout << std::endl;
  algo.echo();
  std::cout << std::endl;
}

/******************************************************************************/
//
// Reads geometry from gmsh grid file
//
/******************************************************************************/
void read_gmsh_grid ( const std::string & grid_fname, Body & aircraft )
{
  unsigned int nnodes, lbl, n, nelements, type, nitems, n0, n1, n2, n3, ntags;
  unsigned int i, j, nthin, nte;
  double x, y, z;
  bool nodes_section, elements_section, thin, te;
  std::string message, line, tag;
  std::ifstream gridfile;
  std::vector<std::string> splitline;

  // Print status

  std::cout << "Reading grid file " + grid_fname + " ..." << std::endl;

  // Open grid file

  gridfile.open(grid_fname.c_str());
  if (not gridfile.is_open())
  {
    message = "Grid file " + grid_fname + " could not be found.";
    conditional_stop(1, "read_gmsh_grid", message);
  }

  // Read until $Nodes header is found

  nodes_section = false;
  while (not nodes_section)
  {
    // Check for end of file

    if (gridfile.eof()) { break;}
    else
    {
      std::getline(gridfile, line);
      if (line.substr(0, 6) == "$Nodes")
      {
        nodes_section = true;
      }
    }
  }

  // Error if required $Nodes section has not been found

  if (not nodes_section)
  {
    conditional_stop(1, "read_gmsh_grid", 
                     "Missing $Nodes section in grid file.");
  }

  // Read number of nodes

  std::getline(gridfile, line);
  nnodes = string2int(line);

  // Read and set nodes in aircraft

  for ( n = 0; n < nnodes; n++ )
  {
    std::getline(gridfile, line);
    splitline = split_string(line);
    lbl = string2int(splitline[0]);
    x = string2double(splitline[1]);
    y = string2double(splitline[2]);
    z = string2double(splitline[3]);
    aircraft.addNode(lbl, x, y, z);
  }

  // Read until $Elements header is found

  elements_section = false;
  while (not elements_section)
  {
    // Check for end of file

    if (gridfile.eof()) { break;}
    else
    {
      std::getline(gridfile, line);
      if (line.substr(0, 9) == "$Elements")
      {
        elements_section = true;
      }
    }
  }

  // Error if required $Elements section has not been found

  if (not elements_section)
  {
    conditional_stop(1, "read_gmsh_grid", 
                     "Missing $Elements section in grid file.");
  }

  // Read number of elements

  std::getline(gridfile, line);
  nelements = string2int(line);

  // Read and set elements in aircraft
  
  nthin = thin_surface_list.size();
  for ( n = 0; n < nelements; n++ )
  {
    std::getline(gridfile, line);
    splitline = split_string(line);
    nitems = splitline.size();
    lbl = string2int(splitline[0]);
    type = string2int(splitline[1]);
    ntags = string2int(splitline[2]);

    // Determine if element is thin surface (no thickness)

    thin = false;
    for ( i = 0; i < ntags; i++ )
    {
      tag = splitline[i+3];  
      for ( j = 0; j < nthin; j++ )
      {
        if (tag == thin_surface_list[j])
        {
          thin = true;
          break;
        }
      }
    }

    // Triangles (type 2)

    if (type == 2)
    {
      n0 = string2int(splitline[nitems-3]);
      n1 = string2int(splitline[nitems-2]);
      n2 = string2int(splitline[nitems-1]);
      if (reverse_normals) { aircraft.addFace(lbl, n2, n1, n0, thin); }
      else { aircraft.addFace(lbl, n0, n1, n2, thin); }
    }

    // Quads (type 3)

    else if (type == 3)
    {
      n0 = string2int(splitline[nitems-4]);
      n1 = string2int(splitline[nitems-3]);
      n2 = string2int(splitline[nitems-2]);
      n3 = string2int(splitline[nitems-1]);
      if (reverse_normals) { aircraft.addFace(lbl, n3, n2, n1, n0, thin); }
      else { aircraft.addFace(lbl, n0, n1, n2, n3, thin); }
    }
  }

  // Close grid file

  gridfile.close();

  // Print status
    
  std::cout << "Successfully read grid file." << std::endl;

  // Set face pointers (generic faces pointing to quadfaces or trifaces)
  
  aircraft.setFacePointers();

  // Construct edges from nodes and faces

  aircraft.constructEdges();

  // Read trailing edges from grid file

  nte = trailing_edge_list.size();
  if (nte > 0)
  {
    // Open grid file
  
    gridfile.open(grid_fname.c_str());
    if (not gridfile.is_open())
    {
      message = "Grid file " + grid_fname + " could not be found.";
      conditional_stop(1, "read_gmsh_grid", message);
    }

    // Read until $Elements header is found
  
    elements_section = false;
    while (not elements_section)
    {
      // Check for end of file
  
      if (gridfile.eof()) { break;}
      else
      {
        std::getline(gridfile, line);
        if (line.substr(0, 9) == "$Elements")
        {
          elements_section = true;
        }
      }
    }
  
    // Error if required $Elements section has not been found
  
    if (not elements_section)
    {
      conditional_stop(1, "read_gmsh_grid", 
                       "Missing $Elements section in grid file.");
    }
  
    // Read number of elements
  
    std::getline(gridfile, line);
    nelements = string2int(line);

    // Read and set trailing edges in aircraft
    
    for ( n = 0; n < nelements; n++ )
    {
      std::getline(gridfile, line);
      splitline = split_string(line);
      nitems = splitline.size();
      lbl = string2int(splitline[0]);
      type = string2int(splitline[1]);
      ntags = string2int(splitline[2]);
  
      // Determine if element is a trailing edge
  
      te = false;
      for ( i = 0; i < ntags; i++ )
      {
        tag = splitline[i+3];  
        for ( j = 0; j < nte; j++ )
        {
          if (tag == trailing_edge_list[j])
          {
            te = true;
            break;
          }
        }
      }
  
      // Lines (type 1)
  
      if (type == 1)
      {
        n0 = string2int(splitline[nitems-2]);
        n1 = string2int(splitline[nitems-1]);
        if (te) { aircraft.setTrailingEdgeFromNodes(n0, n1); }
      }
    }

    // Close grid file

    gridfile.close();
  }
}

/******************************************************************************/
//
// Reads geometry from grid file
//
/******************************************************************************/
void read_grid ( const std::string & grid_format, 
                 const std::string & grid_fname, Body & aircraft )
{
  // Check for recognized file type
  
  if (grid_format == "gmsh") { read_gmsh_grid(grid_fname, aircraft); }
  else
  {
    conditional_stop(1, "read_grid", 
                     "Unrecognized file format " + grid_format + " .");
  }

  // Print info
  
  std::cout << std::endl;
  std::cout << "Body information: " << std::endl;
  std::cout << "Nodes: " << aircraft.numNodes() << std::endl;
  std::cout << "Total edges: " << aircraft.numEdges() << std::endl;
  std::cout << "Trailing edges: " << aircraft.numTrailingEdges() << std::endl;
  std::cout << "Tri faces: " << aircraft.numTris() << std::endl;
  std::cout << "Quad faces: " << aircraft.numQuads() << std::endl;
  std::cout << "Faces belonging to a thin surface: " << aircraft.numThin() 
            << std::endl;
  std::cout << "Faces connected to a trailing edge: " << aircraft.numTEFaces()
            << std::endl;
  if (reverse_normals)
  {
    std::cout << "Note: reversed normals on all faces." << std::endl;
  }
  std::cout << std::endl;
}

/******************************************************************************/
//
// Writes ensight case file
//
/******************************************************************************/
void write_ensight_case ( const std::string & projname,
                          const std::vector<std::string> & outvars,
                          const std::vector<std::string> & outvartypes )
{
  unsigned int space1, space2, nout, i;
  std::string filename, message;
  std::ofstream currfile;

  // Open case file

  filename = projname + ".case";
  currfile.open(filename.c_str());
  if (not currfile.is_open())
  {
    message = "Ensight case file " + filename + " could not be opened.";
    conditional_stop(1, "write_ensight_case", message);
  }

  // Write case file

  space1 = 25;
  space2 = 20;
  nout = outvars.size();
  currfile << "FORMAT" << std::endl;
  currfile << std::setw(space1) << std::left << "type:" 
           << "ensight gold" << std::endl;
  currfile << "GEOMETRY" << std::endl;
  currfile << std::setw(space1) << std::left << "model:" 
           << projname + ".geo" << std::endl;
  currfile << "VARIABLE" << std::endl;
  for ( i = 0; i < nout; i++ )    // Loop over output variables
  {
    if (outvartypes[i] == "elemscalar")
    {
      currfile << std::setw(space1) << std::left << "scalar per element:";
    } 
    else if (outvartypes[i] == "elemvector")
    {
      currfile << std::setw(space1) << std::left << "vector per element:";
    } 
    else if (outvartypes[i] == "nodescalar")
    {
      currfile << std::setw(space1) << std::left << "scalar per node:";
    } 
    else if (outvartypes[i] == "nodevector")
    {
      currfile << std::setw(space1) << std::left << "vector per node:";
    } 
    currfile << std::setw(space2) << std::left << outvars[i];
    currfile << projname + "_" + outvars[i] + ".scl" << std::endl;
  }

  // Close case file

  currfile.close();
}

/******************************************************************************/
//
// Writes ensight geometry file
//
/******************************************************************************/
void write_ensight_geometry ( const std::string & projname, Body & aircraft,
                              Wake & wake, bool edgewrite )
{
  unsigned int tempint, nnodes, ntris, nquads, nedges, nwakelines;
  unsigned int i, j;
  float tempfloat;
  char buffer[80];
  std::string filename, message;
  std::ofstream currfile;

  // Open geometry file

  filename = projname + ".geo";
  currfile.open(filename.c_str(), std::ios::binary);
  if (not currfile.is_open())
  {
    message = "Ensight geometry file " + filename + " could not be opened.";
    conditional_stop(1, "write_ensight_geometry", message);
  }

  // Write geometry file header

  strcpy(buffer, "C Binary");
  currfile.write((char *) & buffer, 80*sizeof(char));
  strcpy(buffer, projname.c_str());
  currfile.write((char *) & buffer, 80*sizeof(char));
  strcpy(buffer, "description line 2");
  currfile.write((char *) & buffer, 80*sizeof(char));
  strcpy(buffer, "node id given");
  currfile.write((char *) & buffer, 80*sizeof(char));
  strcpy(buffer, "element id given");
  currfile.write((char *) & buffer, 80*sizeof(char));

  // Write aircraft part header

  strcpy(buffer, "part");
  currfile.write((char *) & buffer, 80*sizeof(char));
  tempint = 1;                     // Part number
  currfile.write((char *) & tempint, sizeof(unsigned int));
  strcpy(buffer, "Aircraft");
  currfile.write((char *) & buffer, 80*sizeof(char));
  strcpy(buffer, "coordinates");
  currfile.write((char *) & buffer, 80*sizeof(char));

  // Write aircraft nodes

  nnodes = aircraft.numNodes();
  currfile.write((char *) & nnodes, sizeof(unsigned int));
  for ( i = 0; i < nnodes; i++ )   // Write node ids
  {
    tempint = aircraft.node(i).label();
    currfile.write((char *) & tempint, sizeof(unsigned int));
  }
  for ( i = 0; i < nnodes; i++ )   // Write node x-coords
  {
    tempfloat = aircraft.node(i).x();
    currfile.write((char *) & tempfloat, sizeof(float));
  }
  for ( i = 0; i < nnodes; i++ )   // Write node y-coords
  {
    tempfloat = aircraft.node(i).y();
    currfile.write((char *) & tempfloat, sizeof(float));
  }
  for ( i = 0; i < nnodes; i++ )   // Write node z-coords
  {
    tempfloat = aircraft.node(i).z();
    currfile.write((char *) & tempfloat, sizeof(float));
  }

  // Write aircraft tri faces
    
  ntris = aircraft.numTris();
  if (ntris > 0)
  {
    strcpy(buffer, "tria3");
    currfile.write((char *) & buffer, 80*sizeof(char));
    currfile.write((char *) & ntris, sizeof(unsigned int));
    for ( j = aircraft.firstTri(); j <= aircraft.lastTri(); j++ )
    {                              // Write element ids
      tempint = aircraft.face(j).label();
      currfile.write((char *) & tempint, sizeof(unsigned int));
    }
    for ( j = aircraft.firstTri(); j <= aircraft.lastTri(); j++ )
    {                              // Write element nodes
      for ( i = 0; i < 3; i++ )
      {
        tempint = aircraft.face(j).node(i).label();
        currfile.write((char *) & tempint, sizeof(unsigned int));
      }
    }
  }

  // Write aircraft quad faces
    
  nquads = aircraft.numQuads();
  if (nquads > 0)
  {
    strcpy(buffer, "quad4");
    currfile.write((char *) & buffer, 80*sizeof(char));
    currfile.write((char *) & nquads, sizeof(unsigned int));
    for ( j = aircraft.firstQuad(); j <= aircraft.lastQuad(); j++ )
    {                              // Write element ids
      tempint = aircraft.face(j).label();
      currfile.write((char *) & tempint, sizeof(unsigned int));
    }
    for ( j = aircraft.firstQuad(); j <= aircraft.lastQuad(); j++ )
    {                              // Write element nodes
      for ( i = 0; i < 4; i++ )
      {
        tempint = aircraft.face(j).node(i).label();
        currfile.write((char *) & tempint, sizeof(unsigned int));
      }
    }
  }

  // Write wake part header

  strcpy(buffer, "part");
  currfile.write((char *) & buffer, 80*sizeof(char));
  tempint = 2;                     // Part number
  currfile.write((char *) & tempint, sizeof(unsigned int));
  strcpy(buffer, "Wake");
  currfile.write((char *) & buffer, 80*sizeof(char));
  strcpy(buffer, "coordinates");
  currfile.write((char *) & buffer, 80*sizeof(char));

  // Write wake nodes

  nwakelines = wake.numWakeLines();
  nnodes = nwakelines * (num_filaments + 1);
  currfile.write((char *) & nnodes, sizeof(unsigned int));
  tempint = 0;
  for ( i = 0; i < nwakelines; i++ )   // Write node ids
  {
    for ( j = 0; j <= num_filaments; j++ )
    {
      tempint += 1;
      currfile.write((char *) & tempint, sizeof(unsigned int));
    }
  }
  for ( i = 0; i < nwakelines; i++ )   // Write node x-coords
  {
    for ( j = 0; j <= num_filaments; j++ )
    {
      tempfloat = wake.wakeLine(i).x(j);
      currfile.write((char *) & tempfloat, sizeof(float));
    }
  }
  for ( i = 0; i < nwakelines; i++ )   // Write node y-coords
  {
    for ( j = 0; j <= num_filaments; j++ )
    {
      tempfloat = wake.wakeLine(i).y(j);
      currfile.write((char *) & tempfloat, sizeof(float));
    }
  }
  for ( i = 0; i < nwakelines; i++ )   // Write node z-coords
  {
    for ( j = 0; j <= num_filaments; j++ )
    {
      tempfloat = wake.wakeLine(i).z(j);
      currfile.write((char *) & tempfloat, sizeof(float));
    }
  }

  // Write wake edges

  strcpy(buffer, "bar2");
  currfile.write((char *) & buffer, 80*sizeof(char));
  nedges = nwakelines*num_filaments;
  currfile.write((char *) & nedges, sizeof(unsigned int));
  tempint = 0;
  for ( i = 0; i < nwakelines; i++ )   // Write edge ids
  {
    for ( j = 0; j < num_filaments; j++ )
    {
      tempint += 1;
      currfile.write((char *) & tempint, sizeof(unsigned int));
    }
  }
  for ( i = 0; i < nwakelines; i++ )   // Write element nodes
  {
    for ( j = 0; j < num_filaments; j++ )
    {
      tempint = i*(num_filaments+1) + (j + 1);  // First node id for edge
      currfile.write((char *) & tempint, sizeof(unsigned int));
      tempint += 1;                             // Second node id for edge
      currfile.write((char *) & tempint, sizeof(unsigned int));
    }
  }

  // Write edges to new part (if requested)

  if (edgewrite)
  {
    // Header
    
    strcpy(buffer, "part");
    currfile.write((char *) & buffer, 80*sizeof(char));
    tempint = 3;                     // Part number
    currfile.write((char *) & tempint, sizeof(unsigned int));
    strcpy(buffer, "Edges");
    currfile.write((char *) & buffer, 80*sizeof(char));
    strcpy(buffer, "coordinates");
    currfile.write((char *) & buffer, 80*sizeof(char));

    // Nodes
  
    nnodes = aircraft.numNodes();
    currfile.write((char *) & nnodes, sizeof(unsigned int));
    for ( i = 0; i < nnodes; i++ )   // Write node ids
    {
      tempint = aircraft.node(i).label();
      currfile.write((char *) & tempint, sizeof(unsigned int));
    }
    for ( i = 0; i < nnodes; i++ )   // Write node x-coords
    {
      tempfloat = aircraft.node(i).x();
      currfile.write((char *) & tempfloat, sizeof(float));
    }
    for ( i = 0; i < nnodes; i++ )   // Write node y-coords
    {
      tempfloat = aircraft.node(i).y();
      currfile.write((char *) & tempfloat, sizeof(float));
    }
    for ( i = 0; i < nnodes; i++ )   // Write node z-coords
    {
      tempfloat = aircraft.node(i).z();
      currfile.write((char *) & tempfloat, sizeof(float));
    }

    // Edges

    strcpy(buffer, "bar2");
    currfile.write((char *) & buffer, 80*sizeof(char));
    nedges = aircraft.numEdges();
    currfile.write((char *) & nedges, sizeof(unsigned int));
    for ( j = 0; j < nedges; j++ )
    {                              // Write edge ids
      tempint = aircraft.edge(j).index();
      currfile.write((char *) & tempint, sizeof(unsigned int));
    }
    for ( j = 0; j < nedges; j++ )
    {                              // Write element nodes
      for ( i = 0; i < 2; i++ )
      {
        tempint = aircraft.edge(j).node(i).label();
        currfile.write((char *) & tempint, sizeof(unsigned int));
      }
    }
  }

  // Close geometry file

  currfile.close();
}

/******************************************************************************/
//
// Writes ensight element scalar data
//
/******************************************************************************/
void write_ensight_element_scalar ( const std::string & filename,
                                    const std::string & varname, 
                                    Body & aircraft )
{
  unsigned int tempint, ntris, nquads, j;
  float tempfloat;
  char buffer[80];
  std::string message;
  std::ofstream currfile;

  // Open data file

  currfile.open(filename.c_str(), std::ios::binary);
  if (not currfile.is_open())
  {
    message = "Ensight data file " + filename + " could not be opened.";
    conditional_stop(1, "write_ensight_element_scalar", message);
  }

  // Write data file header

  strcpy(buffer, varname.c_str());
  currfile.write((char *) & buffer, 80*sizeof(char));
  strcpy(buffer, "part");
  currfile.write((char *) & buffer, 80*sizeof(char));
  tempint = 1;                     // Part number
  currfile.write((char *) & tempint, sizeof(unsigned int));

  // Write data at tri faces

  ntris = aircraft.numTris();
  if (ntris > 0)
  {
    strcpy(buffer, "tria3");
    currfile.write((char *) & buffer, 80*sizeof(char));
    for ( j = aircraft.firstTri(); j <= aircraft.lastTri(); j++ )
    {
      tempfloat = aircraft.face(j).getScalar(varname);
      currfile.write((char *) & tempfloat, sizeof(float));
    }
  }

  // Write data at quad faces

  nquads = aircraft.numQuads();
  if (nquads > 0)
  {
    strcpy(buffer, "quad4");
    currfile.write((char *) & buffer, 80*sizeof(char));
    for ( j = aircraft.firstQuad(); j <= aircraft.lastQuad(); j++ )
    {
      tempfloat = aircraft.face(j).getScalar(varname);
      currfile.write((char *) & tempfloat, sizeof(float));
    }
  }

  // Close data file

  currfile.close();
}

/******************************************************************************/
//
// Writes ensight element vector data
//
/******************************************************************************/
void write_ensight_element_vector ( const std::string & filename,
                                    const std::string & varname, 
                                    Body & aircraft )
{
  unsigned int tempint, ntris, nquads, j, i;
  char buffer[80];
  float tempfloat;
  std::string message;
  std::ofstream currfile;
  Eigen::Vector3d tempvec;

  // Open data file

  currfile.open(filename.c_str(), std::ios::binary);
  if (not currfile.is_open())
  {
    message = "Ensight data file " + filename + " could not be opened.";
    conditional_stop(1, "write_ensight_element_vector", message);
  }

  // Write data file header

  strcpy(buffer, varname.c_str());
  currfile.write((char *) & buffer, 80*sizeof(char));
  strcpy(buffer, "part");
  currfile.write((char *) & buffer, 80*sizeof(char));
  tempint = 1;                     // Part number
  currfile.write((char *) & tempint, sizeof(unsigned int));

  // Write data at tri faces

  ntris = aircraft.numTris();
  if (ntris > 0)
  {
    strcpy(buffer, "tria3");
    currfile.write((char *) & buffer, 80*sizeof(char));
    for ( i = 0; i < 3; i++ )
    {
      for ( j = aircraft.firstTri(); j <= aircraft.lastTri(); j++ )
      {
        tempfloat = aircraft.face(j).getVector(varname)(i);
        currfile.write((char *) & tempfloat, sizeof(float));
      }
    }
  }

  // Write data at quad faces

  nquads = aircraft.numQuads();
  if (nquads > 0)
  {
    strcpy(buffer, "quad4");
    currfile.write((char *) & buffer, 80*sizeof(char));
    for ( i = 0; i < 3; i++ )
    {
      for ( j = aircraft.firstQuad(); j <= aircraft.lastQuad(); j++ )
      {
        tempfloat = aircraft.face(j).getVector(varname)(i);
        currfile.write((char *) & tempfloat, sizeof(float));
      }
    }
  }

  // Close data file

  currfile.close();
}

/******************************************************************************/
//
// Writes ensight node scalar data
//
/******************************************************************************/
void write_ensight_node_scalar ( const std::string & filename,
                                 const std::string & varname, 
                                 Body & aircraft )
{
}

/******************************************************************************/
//
// Writes ensight node vector data
//
/******************************************************************************/
void write_ensight_node_vector ( const std::string & filename,
                                 const std::string & varname, 
                                 Body & aircraft )
{
}

/******************************************************************************/
//
// Writes solution output files in ensight format
//
/******************************************************************************/
void write_ensight_output ( const std::string & projname,
                            const std::vector<std::string> & outvars,
                            const std::vector<std::string> & outvartypes,
                            Body & aircraft, Wake & wake )
{
  unsigned int n, nout;
  std::string filename;

  // Write case file

  write_ensight_case(projname, outvars, outvartypes);

  // Write geometry file

  write_ensight_geometry(projname, aircraft, wake, write_edges);

  // Write variable files

  nout = outvars.size();
  for ( n = 0; n < nout; n++ )
  {
    filename = projname + "_" + outvars[n] + ".scl";
    if (outvartypes[n] == "elemscalar")
    {
      write_ensight_element_scalar(filename, outvars[n], aircraft);
    }
    else if (outvartypes[n] == "elemvector")
    {
      write_ensight_element_vector(filename, outvars[n], aircraft);
    }
    else if (outvartypes[n] == "nodescalar")
    {
      write_ensight_node_scalar(filename, outvars[n], aircraft);
    }
    else if (outvartypes[n] == "nodevector")
    {
      write_ensight_node_vector(filename, outvars[n], aircraft);
    }
  }
}

/******************************************************************************/
//
// Writes solution output files
//
/******************************************************************************/
void write_solution_output ( const std::string & format,
                             const std::string & projname,
                             const std::vector<std::string> & outvars,
                             Body & aircraft, Wake & wake )
{
  unsigned int i, j, nout, nsupported;
  std::vector<std::string> supportedvars, supportedvartypes, outvartypes;
  bool found;

  // Print notification

  std::cout << "Writing solution data to " + format + " format ..." 
            << std::endl;

  // Create index vectors of supported output variables and types

  supportedvars.push_back("area");
  supportedvartypes.push_back("elemscalar");
  supportedvars.push_back("normal");
  supportedvartypes.push_back("elemvector");

  nsupported = supportedvars.size();

  // Associate requested output variables with correct type

  nout = outvars.size();
  outvartypes.resize(nout);
  for ( i = 0; i < nout; i++ )
  {
    found = false;
    for ( j = 0; j < nsupported; j++ )
    {
      if (outvars[i] == supportedvars[j])
      {
        found = true;
        outvartypes[i] = supportedvartypes[j];
        break;
      }
    }
    if (not found)
    {
      conditional_stop(1, "write_solution_output",
            "Requested output variable '" + outvars[i] + "' is not supported.");
    }
  } 

  // Check for recognized file type and write solution files

  if (format == "ensight") 
  { 
    write_ensight_output(projname, outvars, outvartypes, aircraft, wake); 
  }
  else
  {
    conditional_stop(1, "write_solution_output",
                     "Unrecognized file format " + format + ".");
  }

  // Print notification

  std::cout << "Finished writing solution data." << std::endl;
  std::cout << std::endl;
}
