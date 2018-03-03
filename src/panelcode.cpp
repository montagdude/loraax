#include <fenv.h>
#include "util.h"
#include "settings.h"
#include "body.h"
#include "wake.h"
#include "wake_relaxation.h"
#include "file_io.h"

int main (int argc, char* argv[]) 
{ 
  Body aircraft;
  Wake wake;

  // Get input file name from command line arguments
  
  if (argc < 2)
  {
    conditional_stop(1, "main", 
                     "Input file name required as command line argument.");
  }

  // Catch floating point exceptions

  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  // Read inputs from file

  read_settings(argv[1]);

  // Read geometry from grid file

  read_grid(grid_format, grid_file, aircraft);

  // Initialize wake

  initialize_wake(aircraft, wake);

  // Write output data to file

  write_solution_output(output_file_format, case_name,
                        output_variable_list, aircraft, wake);

  return 0; 
}
