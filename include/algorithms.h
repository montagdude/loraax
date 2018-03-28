// Contains various algorithms

#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <vector>

class SectionalObject;

// Search settings

struct simplex_options_type
{
  double tol;
  unsigned int maxit;
  bool display_progress;
};

struct conjgrad_options_type
{
  double tol, h, dxmax;
  unsigned int maxit;
  bool display_progress;
}; 

// Public routines

std::vector<double> uniform_spacing ( const double & slen, unsigned int n );
std::vector<double> cosine_spacing ( const double & slen, unsigned int n );
std::vector<double> sine_spacing ( const double & slen, unsigned int n );
void opt_tanh_spacing ( unsigned int n, const double & slen, const double & sp0,
                        const double & sp1, double & a4, double & a5 );
double tanh_spacing ( unsigned int i, const double & a4, const double & a5,
                      unsigned int n, const double & slen, const double & sp0,
                      const double & sp1 );
void simplex_search ( std::vector<double> & xopt, double & fmin,
                      unsigned int & nsteps, unsigned int & fevals,
                      double (*objfunc)(const std::vector<double> & x),
                      const std::vector<double> & x0,
                      const simplex_options_type & searchopt );
void sort_sections ( SectionalObject *sections[], unsigned int nsecs );
void golden_search ( double & xmin, double & fmin, const double bounds[2],
                     unsigned int & fevals, const double & tol,
                     const unsigned int itmax,
                     double (*objfunc)(const double & x) );
void gradient ( const std::vector<double> & x, std::vector<double> & grad,
                double (*func)(const std::vector<double> & x),
                const double & h=1.E-08 );
void conjgrad_search ( std::vector<double> & xopt, double & fmin,
                       unsigned int & nsteps, unsigned int & fevals,
                       double (*objfunc)(const std::vector<double> & x),
                       const std::vector<double> & x0,
                       const conjgrad_options_type & searchopt );

#endif
