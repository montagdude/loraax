#include <iostream>
#include <iomanip>
#include <string>
#include <ctime>
#include <cmath>
#include <sys/stat.h>
#include <dirent.h>
#include "clo_parser.h"
#include "settings.h"
#include "aircraft.h"
#include "util.h"

#ifndef LORAAX_VERSION
    #define LORAAX_VERSION ""
#endif

void create_or_backup_dir ( const std::string & dirname )
{
    DIR *pdir = NULL;
    time_t now;
    tm *ltm;
    std::string newpath;
    
    pdir = opendir(dirname.c_str());
    if (pdir != NULL)
    {
        pdir = opendir("backup");
        if (pdir == NULL)
#ifdef ISMINGW
            mkdir("backup");
#else
            mkdir("backup", S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
#endif
        now = time(0);
        ltm = localtime(&now);
        newpath = "backup/" + dirname;
        newpath += "." + int2string(1900 + ltm->tm_year);
        newpath += "." + int2string(1 + ltm->tm_mon);
        newpath += "." + int2string(ltm->tm_mday);
        newpath += "." + int2string(ltm->tm_hour);
        newpath += "." + int2string(ltm->tm_min);
        newpath += "." + int2string(ltm->tm_sec);
        rename(dirname.c_str(), newpath.c_str());
    }
    
#ifdef ISMINGW
    mkdir(dirname.c_str());
#else
    mkdir(dirname.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
#endif
}

int main (int argc, char* argv[])
{
    CLOParser parser;
    std::string geom_file;
    int check;
    Aircraft ac;
    unsigned int iter, viz_iter;
    double lift, oldlift;
    bool converged;
    
    // Parse CLOs
    
    check = parser.checkCLOs(argc, argv, LORAAX_VERSION);
    if (check < 0)
      return 0;
    else if (check > 0)
      return 1;
    else
    {
      // Read settings
    
        std::cout << "Reading settings ..." << std::endl;
        if (read_settings(parser.inputFile(), geom_file) != 0)
            return 2;
    }
    std::cout << "Freestream Mach: "
              << std::setprecision(5) << std::setw(8) << std::left << minf
              << std::endl;
    
    // Read geometry
    
    std::cout << "Reading and discretizing geometry ..." << std::endl;
    if (ac.readXML(geom_file) != 0)
        return 3;
    
    // Set up output directories or back up existing
    
    create_or_backup_dir("visualization");
    create_or_backup_dir("sectional");
    create_or_backup_dir("forcemoment");
    
    // Iterate
    
    iter = 0;
    viz_iter = 0;
    converged = false;
    while (int(iter) < maxiters)
    {
        iter++;
        viz_iter++;

      // Convect wake
    
        std::cout << "Iteration " << iter << std::endl;
        
        if ( (iter > 1) && rollup_wake )
        {
            std::cout << "  Convecting wake ..." << std::endl;
            ac.moveWake();
        }
        
        // Set source strengths
        
        std::cout << "  Setting source strengths ..." << std::endl;
        ac.setSourceStrengths(iter==1);
        
        // Construct, factorize, and solve the system
        
        std::cout << "  Constructing the linear system ..." << std::endl;
        ac.constructSystem(iter==1);
        if ( (iter == 1) || rollup_wake )
        {
            std::cout << "  Factorizing the AIC matrix ..." << std::endl;
            ac.factorize();
        }
        std::cout << "  Solving the linear system with " << ac.systemSize()
                  << " unknowns ..." << std::endl;
        ac.solveSystem();
        
        // Set doublet strengths on surface and wake
        
        std::cout << "  Setting doublet strengths ..." << std::endl;
        ac.setDoubletStrengths();
        
        // Compute surface velocities and pressures
        
        std::cout << "  Computing surface velocity and pressure ..."
                  << std::endl;
        ac.computeSurfaceQuantities();
        
        // Viscous BL computations with Xfoil
        
        if (viscous)
        {
            std::cout << "  Computing viscous BL with Xfoil ..." << std::endl;
            ac.computeBL();
            if (iter == 1)
                ac.setupViscousWake();
        }
        
        // Compute forces and moments
        
        std::cout << "  Computing forces and moments ..." << std::endl;
        ac.computeForceMoment();
        ac.writeForceMoment(iter);
        lift = ac.liftCoefficient();
        if (viscous)
        {
            std::cout << "  CL: " << std::setprecision(5) << std::setw(8)
                      << std::left << ac.liftCoefficient();
            std::cout << "  CD: " << std::setprecision(5) << std::setw(8)
                      << std::left << ac.dragCoefficient() << std::endl;
            std::cout << "  CDinduced: " << std::setprecision(5)
                      << std::setw(8) << std::left
                      << ac.inducedDragCoefficient();
            std::cout << "  CDparasitic: " << std::setprecision(5)
                      << std::setw(8) << std::left
                      << ac.parasiticDragCoefficient() << std::endl;
            std::cout << "  Cm: " << std::setprecision(5) << std::setw(8)
                      << std::left << ac.pitchingMomentCoefficient()
                      << std::endl;
        }
        else
        {
            std::cout << "  CL: " << std::setprecision(5) << std::setw(8)
                      << std::left << ac.trefftzLiftCoefficient();
            std::cout << "  CD: " << std::setprecision(5) << std::setw(8)
                      << std::left << ac.inducedDragCoefficient();
            std::cout << "  (Trefftz plane)" << std::endl;
            std::cout << "  CL: " << std::setprecision(5) << std::setw(8)
                      << std::left << ac.integratedLiftCoefficient();
            std::cout << "  CD: " << std::setprecision(5) << std::setw(8)
                      << std::left << ac.integratedDragCoefficient();
            std::cout << "  (Pressure integration)" << std::endl;
            std::cout << "  Cm: " << std::setprecision(5) << std::setw(8)
                      << std::left << ac.pitchingMomentCoefficient()
                      << std::endl;
        }
        
        // Write visualization
        
        if ( (int(viz_iter) == viz_freq) && (int(iter) < maxiters) )
        {
            std::cout << "  Writing VTK visualization ..." << std::endl;
            ac.writeViz(casename, iter);
            ac.writeSectionForceMoment(iter);
            viz_iter = 0;
        }

        // Stop iterating if converged

        if ( (iter > 1) && (int(iter) >= miniters) )
        {
            if (std::abs(lift - oldlift) < stop_tol)
            {
                std::cout << "Solution is converged." << std::endl;
                converged = true;
                break;
            }
        }
        else if ( (! rollup_wake) && (! viscous) )
        {
            std::cout << "Solution is converged." << std::endl;
            converged = true;
            break;
        }

        oldlift = lift;
    }

    if (! converged)
        print_warning("main", "Solution did not converge.");
    
    // Write final visualization
    
    if (int(viz_iter) != viz_freq)
    {
        std::cout << "Writing final VTK visualization ..." << std::endl;
        ac.writeViz(casename, iter);
        ac.writeSectionForceMoment(iter);
    }
    
    return 0;
}
