// Contains various algorithms

#ifndef ALGORITHMS_H
#define ALGORITHMS_H

#include <vector>

// Simplex search settings

struct simplex_options_type
{
  double tol;
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

#endif
