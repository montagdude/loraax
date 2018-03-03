// Stores and reads global settings

#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>
#include <vector>
#include <Eigen/Core>

// Global variables (settings)

extern std::string grid_format;          // Grid file format
extern std::string grid_file;            // Grid filename
extern bool reverse_normals;             // Reverses all face normals on read

extern double alpha;                     // Angle of attack (degrees)
extern double beta;                      // Angle of sideslip (degrees)
extern Eigen::Vector3d uinf;             // Freestream velocity (computed)

extern double ref_area;                  // Force and moment reference area
extern double mom_len_x, mom_len_y, mom_len_z;
                                         // Reference lengths for moment 
                                         //   coefficients
extern double mom_cen_x, mom_cen_y, mom_cen_z;
                                         // Location of moment calculation
                                   
extern std::vector<std::string> trailing_edge_list;   
                                         // List of trailing edges
extern std::vector<std::string> thin_surface_list;   
                                         // List of thin surfaces

extern double deforming_wake_length;     // Length of deforming/visualized wake
                                         //   (last filament goes to infinity)
extern unsigned int num_filaments;       // Number of filaments in each line
extern double core_radius_factor;        // Vortex filament core radius/length

extern std::string case_name;            // Case name
extern std::string output_file_format;   // Output file format
extern std::vector<std::string> output_variable_list;
                                         // List of output variables
extern bool write_edges;                 // Whether to write edges in viz.

extern double farfield_distance_factor;  // Distance from panel at which point
                                         //   singularities are used in place of
                                         //   distributed singularities
extern bool relax_wake;                  // Whether to relax wake
extern double wake_convergence_tolerance;// RMS of wake node positions before
                                         //   it is considered converged

// Public functions
  
void read_settings ( const std::string & ); 
                                         // Reads settings from namelist file

#endif
