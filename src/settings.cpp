// Stores global settings

#include <string>
#include <vector>
#include <Eigen/Core>

std::string grid_format;           // Grid file format
std::string grid_file;             // Grid filename
bool reverse_normals;              // Reverses all face normals on read

double alpha;                      // Angle of attack (degrees)
double beta;                       // Angle of sideslip (degrees)
Eigen::Vector3d uinf;              // Freestream velocity (computed)

double ref_area;                   // Force and moment reference area
double mom_len_x, mom_len_y, mom_len_z;
                                   // Reference lengths for moment coefficients
double mom_cen_x, mom_cen_y, mom_cen_z;
                                   // Location of moment calculation
                                     
std::vector<std::string> trailing_edge_list;    
                                   // List of trailing edges
std::vector<std::string> thin_surface_list;    
                                   // List of thin surfaces

double deforming_wake_length;      // Length of deforming / visualized wake
                                   //   (last filament goes to infinity)
unsigned int num_filaments;        // Number of filaments in each line
double core_radius_factor;         // Vortex core radius / min te width

std::string case_name;             // Case name
std::string output_file_format;    // Output file format
std::vector<std::string> output_variable_list;
                                   // List of output variables
bool write_edges;                  // Whether to write edges in visualization

double farfield_distance_factor;   // Distance from panel at which point
                                   //   singularities are used in place of
                                   //   distributed singularities
bool relax_wake;                   // Whether to relax wake
double wake_convergence_tolerance; // RMS of wake node positions before it is
                                   //   considered converged
