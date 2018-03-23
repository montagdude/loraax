// Contains various algorithms

#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>	// min_element
#include "util.h"
#include "sectional_object.h"
#include "algorithms.h"

// Global parameters for tanh spacing objective function

unsigned int tanh_N;
double tanh_slen, tanh_sp0, tanh_sp1;

/******************************************************************************/
//
// Uniform spacing function
//
/******************************************************************************/
std::vector<double> uniform_spacing ( const double & slen, unsigned int n )
{
  unsigned int i;
  std::vector<double> svec(n, 0.0);

  for ( i = 0; i < n; i++ )
  {
    svec[i] = double(i) / double(n-1) * slen;
  }

  return svec; 
}

/******************************************************************************/
//
// Cosine spacing function: clustering at x = 0 and x = slen
//
/******************************************************************************/
std::vector<double> cosine_spacing ( const double & slen, unsigned int n )
{
  double theta;
  unsigned int i;
  std::vector<double> svec(n, 0.0);

  for ( i = 0; i < n; i++ )
  {
    theta = (double(i) / double(n-1) - 1.) * M_PI;
    svec[i] = 0.5*(cos(theta) + 1.0) * slen;
  }

  return svec; 
}

/******************************************************************************/
//
// Sine spacing function: clustering at x = slen
//
/******************************************************************************/
std::vector<double> sine_spacing ( const double & slen, unsigned int n )
{
  double theta;
  unsigned int i;
  std::vector<double> svec(n, 0.0);

  for ( i = 0; i < n; i++ )
  {
    theta = double(i) / double(n-1) * 0.5 * M_PI;
    svec[i] = sin(theta) * slen;
  }

  return svec; 
}

/******************************************************************************/
//
// tanh spacing basis functions corresponding to a4 and a5 coefficients
//
/******************************************************************************/
double phi4 ( unsigned int i, unsigned int N )
{
  return tanh(2.*double(i)/double(N));
}

double phi5 ( unsigned int i, unsigned int N )
{
  return tanh(2.*(1. - double(i)/double(N)));
}

/******************************************************************************/
//
// Computes stretching for tanh spacing function. Adds penatly for any panels
// with length <= 0.
//
/******************************************************************************/
double tanh_stretching ( const std::vector<double> & coef )
{
  unsigned int N, i;
  double sp0, sp1, slen;
  double a1, a2, a3, a4, a5;
  double a1sum, a2sum, a3sum, a4sum, a5sum;
  double dbleN, maxstretch, stretch, penaltyval, dm, dp;
  std::vector<double> rhs(2, 0.);

  // Set coefficients and other parameters

  a4 = coef[0];
  a5 = coef[1];
  N = tanh_N;
  sp0 = tanh_sp0;
  sp1 = tanh_sp1;
  slen = tanh_slen;

  // Solve for a1, a2, and a3 to satisfy constraints on slen, sp0, and sp1

  a1 = sp0 - a4*phi4(0,N) - a5*phi5(0,N);
  a1sum = 0.;
  a2sum = 0.;
  a3sum = 0.;
  a4sum = 0.;
  a5sum = 0.;
  for ( i = 0; i <= N; i++ )
  {
    a1sum += 1.;
    a2sum += double(i);
    a3sum += pow(double(i),2.);
    a4sum += phi4(i,N);
    a5sum += phi5(i,N);
  }

  rhs[0] = sp1 - a1 - a4*phi4(N,N) - a5*phi5(N,N);
  rhs[1] = slen - a1*a1sum - a4*a4sum - a5*a5sum;
  dbleN = double(N);
  a3 = (rhs[1] - a2sum*rhs[0]/dbleN ) / (a3sum - a2sum*dbleN);
  a2 = (rhs[0] - a3*pow(dbleN,2.)) / dbleN;

  // Compute max stretching and add penalty for any panels with length <= 0

  maxstretch = 0.;
  penaltyval = 0.;
  for ( i = 1; i <= N; i++ )
  {
    dm = a1 + a2*double(i-1) + a3*pow(double(i-1),2.) + a4*phi4(i-1,N)
       + a5*phi5(i-1,N);
    dp = a1 + a2*double(i) + a3*pow(double(i),2.) + a4*phi4(i,N) + a5*phi5(i,N);
    if (dm > dp)
      stretch = (dm - dp)/dp;
    else
      stretch = (dp - dm)/dm;
    if (stretch > maxstretch)
      maxstretch = stretch;

    if (dm <= 0.)
      penaltyval += 1.;
    if ( (i == N) && (dp <= 0.) )
      penaltyval += 1.;
  }

  return maxstretch + penaltyval;
}

/******************************************************************************/
//
// Computes tanh spacing coefficients a4 and a5 to minimize stretchting. See
// tanh_spacing for description of inputs and outputs.
//
/******************************************************************************/
void opt_tanh_spacing ( unsigned int n, const double & slen, const double & sp0,
                        const double & sp1, double & a4, double & a5 )
{
  std::vector<double> optcoef, coef0;
  double maxstretch;
  unsigned int steps, fevals;
  simplex_options_type searchopt;
  double (*objfunc)(const std::vector<double> & x);

  // Simplex search options

  searchopt.tol = 1.E-12;
  searchopt.maxit = 1000;
#ifdef DEBUG
  searchopt.display_progress = true;
#else
  searchopt.display_progress = false;
#endif

  // Global parameters for tanh spacing

  tanh_N = n - 2;
  tanh_slen = slen;
  tanh_sp0 = sp0;
  tanh_sp1 = sp1;

  // Compute optimal tanh spacing coefficients

  optcoef.resize(2);
  coef0.resize(2);
  coef0[0] = 0.;
  coef0[1] = 0.;
  objfunc = &tanh_stretching;
  simplex_search(optcoef, maxstretch, steps, fevals, objfunc, coef0, searchopt);
  a4 = optcoef[0];
  a5 = optcoef[1];
}  

/******************************************************************************/
//
// Computes local point spacing using tanh spacing function:
// space[i] = a1 + a2*i + a3*i^2 + a4*tanh(2*i/(N-2)) + a5*tanh(2*(1-i/(N-2)))
// 
// Inputs:
//   i: current index [0, n-2]
//   a4, a5: tanh coefficients (a1 - a3 calculated automatically)
//   n: number of points on curve, equal to N+2
//   slen: curve length
//   s0, s1: initial and final spacings
//
/******************************************************************************/
double tanh_spacing ( unsigned int i, const double & a4, const double & a5,
                      unsigned int n, const double & slen, const double & sp0,
                      const double & sp1 )
{
  unsigned int N, j;
  double a1, a2, a3;
  double a1sum, a2sum, a3sum, a4sum, a5sum;
  double dbleN;
  double rhs[2];

  N = n - 2;
  
  // Solve for a1, a2, and a3 to satisfy constraints on slen, sp0, and sp1

  a1 = sp0 - a4*phi4(0,N) - a5*phi5(0,N);
  a1sum = 0.;
  a2sum = 0.;
  a3sum = 0.;
  a4sum = 0.;
  a5sum = 0.;
  for ( j = 0; j <= N; j++ )
  {
    a1sum += 1.;
    a2sum += double(j);
    a3sum += pow(double(j),2.);
    a4sum += phi4(j,N);
    a5sum += phi5(j,N);
  }

  rhs[0] = sp1 - a1 - a4*phi4(N,N) - a5*phi5(N,N);
  rhs[1] = slen - a1*a1sum - a4*a4sum - a5*a5sum;
  dbleN = double(N);
  a3 = (rhs[1] - a2sum*rhs[0]/dbleN ) / (a3sum - a2sum*dbleN);
  a2 = (rhs[0] - a3*pow(dbleN,2.)) / dbleN;

  return a1 + a2*double(i) + a3*pow(double(i),2.) + a4*phi4(i,N) + a5*phi5(i,N);
}

/******************************************************************************/
//
// Sorts designs according to their objective function values
//
/******************************************************************************/
void sort_designs ( std::vector<std::vector<double> > & dv,
                    std::vector<double> & objvals )
{
  std::vector<std::vector<double> > tempdv;
  std::vector<double> tempvals;
  std::vector<unsigned int> finalorder, temporder;
  unsigned int i, ndes, nvars, sortcounter;
  bool sorted;

  // Set vector sizes

  ndes = objvals.size();
  nvars = dv[0].size();
  tempdv.resize(ndes);
  for ( i = 0; i < ndes; i++ )
  {
    tempdv[i].resize(nvars);
  }
  tempvals.resize(ndes);
  finalorder.resize(ndes);
  temporder.resize(ndes);

#ifdef DEBUG
  if (dv.size() != ndes)
    conditional_stop(1, "sort_designs", "dv and objvals size mismatch.");
#endif

  // Set up indexing array

  for ( i = 0; i < ndes; i++ )
  {
    finalorder[i] = i;
  }
  temporder = finalorder;

  // Bubble sorting algorithm

  sorted = false;
  tempvals = objvals;
  while (! sorted)
  {
    sortcounter = 0;
    for ( i = 0; i < ndes-1; i++ )
    {
      if (objvals[i+1] < objvals[i])
      {
        // Flip order

        tempvals[i] = objvals[i+1];
        tempvals[i+1] = objvals[i];
        temporder[i] = finalorder[i+1];
        temporder[i+1] = finalorder[i];

        finalorder[i] = temporder[i];
        finalorder[i+1] = temporder[i+1];
        objvals[i] = tempvals[i];
        objvals[i+1] = tempvals[i+1];
        sortcounter += 1;
      }
    }
    if (sortcounter == 0)
      sorted = true;
  }

  // Use indexing array to rearrange designs

  for ( i = 0; i < ndes; i++ )
  {
    tempdv[i] = dv[finalorder[i]];
  }
  dv = tempdv;
}

/******************************************************************************/
//
// Nelder-Mead simplex search algorithm
//
/******************************************************************************/
void simplex_search ( std::vector<double> & xopt, double & fmin,
                      unsigned int & nsteps, unsigned int & fevals,
                      double (*objfunc)(const std::vector<double> & x),
                      const std::vector<double> & x0,
                      const simplex_options_type & searchopt )
{
  std::vector<std::vector<double> > dv;
  std::vector<double> objvals, xcen, xr, xe, xc;
  double rho, xi, gam, sigma, fr, fe, fc, dist, diam;
  unsigned int i, j, k, nvars;
  bool converged, needshrink;

  // Standard Nelder-Mead constants

  rho = 1.;
  xi = 2.;
  gam = 0.5;
  sigma = 0.5;

  // Set vector sizes

  nvars = x0.size();
  dv.resize(nvars+1);
  for ( i = 0; i < nvars+1; i++ )
  {
    dv[i].resize(nvars);
  }
  objvals.resize(nvars+1);
  xcen.resize(nvars);
  xr.resize(nvars);
  xe.resize(nvars);
  xc.resize(nvars);

  // Set up initial simplex

  if (searchopt.display_progress)
    std::cout << "Initializing simplex ..." << std::endl;
  fevals = 0;
  for ( j = 0; j < nvars; j++ )
  {
    for ( i = 0; i < nvars; i++ )
    {
      if (i == j)
      {
        if (x0[i] == 0.)
          dv[j][i] = 0.00025;
        else
          dv[j][i] = 1.05*x0[i];
      }
      else
        dv[j][i] = x0[i];
    }
    objvals[j] = objfunc(dv[j]);
    fevals += 1;
  } 
  dv[nvars] = x0;
  objvals[nvars] = objfunc(x0);
  fevals += 1;
  fmin = vector_min(objvals);

  nsteps = 0;
  needshrink = false;
  converged = false;
  if (searchopt.display_progress)
    std::cout << "Simplex optimization progress:" << std::endl;
  while (! converged)
  {
    nsteps += 1;
    if (nsteps == searchopt.maxit)
      converged = true;

    // Sort according to ascending objective function value

    sort_designs(dv, objvals);
    fmin = objvals[0];

    // Compute max distance between designs

    diam = 0;
    for ( j = 1; j < nvars+1; j++ )
    {
      dist = 0.;
      for ( k = 0; k < nvars; k++ )
      {
        dist += pow(dv[0][k]-dv[j][k], 2.);
      }
      dist = sqrt(dist);
      if (dist > diam)
        diam = dist;
    }

    // Check for convergence and display progress

    if (diam < searchopt.tol)
      converged = true;

    if (searchopt.display_progress)
      std::cout << "  Iteration: " << nsteps << " min objfunc value: " << fmin
                << std::endl;

    // Compute centroid of best nvars designs

    std::fill(xcen.begin(), xcen.end(), 0.);
    for ( j = 0; j < nvars; j++ )
    {
      for ( i = 0; i < nvars; i++ )
      {
        xcen[i] += dv[j][i];
      }
    }
    for ( i = 0; i < nvars; i++ )
    {
      xcen[i] /= double(nvars);
    }

    // Compute the reflection point and evaluate its objective function value

    for ( i = 0; i < nvars; i++ )
    {
      xr[i] = (1.+rho)*xcen[i] - rho*dv[nvars][i];
    }
    fr = objfunc(xr);
    fevals += 1;

    if ( (objvals[0] <= fr) && (fr < objvals[nvars-1]) )
    {
      // Accept reflection point

      dv[nvars] = xr;
      objvals[nvars] = fr;
      continue;
    }
    else if (fr < objvals[0])
    {
      // Expand

      for ( i = 0; i < nvars; i++ )
      {
        xe[i] = (1.+rho*xi)*xcen[i] - rho*xi*dv[nvars][i];
      }
      fe = objfunc(xe);
      fevals += 1;
      if (fe < fr)
      {
        dv[nvars] = xe;
        objvals[nvars] = fe;
      }
      else
      {
        dv[nvars] = xr;
        objvals[nvars] = fr;
      }
      continue;
    }
    else if (fr >= objvals[nvars-1])
    {
      if (fr < objvals[nvars])
      {
        // Outside contraction

        for ( i = 0; i < nvars; i++ )
        {
          xc[i] = (1.+rho*gam)*xcen[i] - rho*gam*dv[nvars][i];
        }
        fc = objfunc(xc);
        fevals += 1;

        if (fc < objvals[nvars])
        {
          dv[nvars] = xc;
          objvals[nvars] = fc;
          needshrink = false;
        }
        else
          needshrink = true;
      }
      else
      {
        // Inside contraction

        for ( i = 0; i < nvars; i++ )
        {
          xc[i] = (1.-gam)*xcen[i] + gam*dv[nvars][i];
        }
        fc = objfunc(xc);
        fevals += 1;

        if (fc < objvals[nvars])
        {
          dv[nvars] = xc;
          objvals[nvars] = fc;
          needshrink = false;
        }
        else
          needshrink = true;
      }

      if (needshrink)
      {
        // Shrink

        for ( j = 1; j < nvars+1; j++ )
        {
          for ( i = 0; i < nvars; i++ )
          {
            dv[j][i] = dv[j][0] + sigma*(dv[j][i] - dv[j][0]);
          }
          objvals[j] = objfunc(dv[j]);
          fevals += 1;
        }
        continue;
      }
    }     
  }

  // Sort final result according to ascending objective function value

  sort_designs(dv, objvals);
  xopt = dv[0];
  fmin = objvals[0];

  // Compute max distance between designs

  diam = 0;
  for ( j = 1; j < nvars+1; j++ )
  {
    dist = 0.;
    for ( k = 0; k < nvars; k++ )
    {
      dist += pow(dv[0][k]-dv[j][k], 2.);
    }
    dist = sqrt(dist);
    if (dist > diam)
      diam = dist;
  }

  // Display warning if max iterations are reached

  if ( (nsteps == searchopt.maxit) && (diam >= searchopt.tol) )
  {
    print_warning("simplex_search", "Max number of iterations was reached.");
  }
}

/******************************************************************************/
//
// Sorts sections based on y
//
/******************************************************************************/
void sort_sections ( SectionalObject *sections[], unsigned int nsecs )
{
  SectionalObject *tempsections[nsecs];
  unsigned int i, sortcounter;
  bool sorted = false;

  // Bubble sorting algorithm

  sorted = false;
  for ( i = 0; i < nsecs; i++ )
  {
    tempsections[i] = sections[i];
  }
  while (! sorted)
  {
    sortcounter = 0;
    for ( i = 0; i < nsecs-1; i++ )
    {
      if (sections[i+1]->y() < sections[i]->y())
      {
        // Flip order

        tempsections[i] = sections[i+1];
        tempsections[i+1] = sections[i];
        sections[i] = tempsections[i];
        sections[i+1] = tempsections[i+1];
        sortcounter += 1;
      }
    }
    if (sortcounter == 0)
      sorted = true;
  }
}
