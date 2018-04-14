// Computes velocity potential and velocity due to singularity elements

#define _USE_MATH_DEFINES

#include <cmath>
#include <algorithm>	// std::max
#include <string>
#include <Eigen/Core>
#include "util.h"

// Note: all distributed source/doublet distributions on panels assume that 
// panel centroid is at the origin and that the panel endpoints lie on the z = 0
// plane.

const double eps = 1.E-12;       
const double corefact = 0.025;   // Scaling factor for vortex core radius.
                                 //   0.025 results in error between 1/d and f1 
                                 //   less than 5% when d = rcore.

/******************************************************************************/
//
// Computes potential at a point due to a point source with unit strength. 
// Source is assumed to be at the origin.
//
/******************************************************************************/
double point_source_potential ( const double & x, const double & y,
                                const double & z )
{
  return -1. / (4.*M_PI*std::sqrt(x*x + y*y + z*z));
}

/******************************************************************************/
//
// Computes potential at a point due to triangular source panel with
// unit strength. Panel endpoints given in clockwise order.
//
/******************************************************************************/
double tri_source_potential ( 
                        const double & x, const double & y, const double & z,
                        const double & x0, const double & y0, const double & x1,
                        const double & y1, const double & x2, const double & y2,
                        bool onpanel, const std::string & side )
{
  double d01, d12, d20;
  double r0, r1, r2;
  double dx01, dx12, dx20, m01, m12, m20;
  double e0, e1, e2;
  double h0, h1, h2;
  double phi, myz;

  if (onpanel)
    myz = 0.;
  else
    myz = z;

  // Geometric quantities needed

  d01 = std::sqrt(std::pow(x1 - x0, 2.) + std::pow(y1 - y0, 2.));
  d12 = std::sqrt(std::pow(x2 - x1, 2.) + std::pow(y2 - y1, 2.));
  d20 = std::sqrt(std::pow(x0 - x2, 2.) + std::pow(y0 - y2, 2.));

  r0 = std::sqrt(std::pow(x - x0, 2.) + std::pow(y - y0, 2.) + 
                 std::pow(myz, 2.));
  r1 = std::sqrt(std::pow(x - x1, 2.) + std::pow(y - y1, 2.) + 
                 std::pow(myz, 2.));
  r2 = std::sqrt(std::pow(x - x2, 2.) + std::pow(y - y2, 2.) + 
                 std::pow(myz, 2.));

  dx01 = sign(x1-x0)*std::max(std::abs(x1-x0), eps);
  dx12 = sign(x2-x1)*std::max(std::abs(x2-x1), eps);
  dx20 = sign(x0-x2)*std::max(std::abs(x0-x2), eps);
  m01 = (y1 - y0)/dx01;
  m12 = (y2 - y1)/dx12;
  m20 = (y0 - y2)/dx20;

  e0 = std::pow(x - x0, 2.) + std::pow(myz, 2.);
  e1 = std::pow(x - x1, 2.) + std::pow(myz, 2.);
  e2 = std::pow(x - x2, 2.) + std::pow(myz, 2.);

  h0 = (x - x0)*(y - y0);
  h1 = (x - x1)*(y - y1);
  h2 = (x - x2)*(y - y2);

  // Source potential

  phi = -1. / (4.*M_PI) * (
          ((x-x0)*(y1-y0) - (y-y0)*(x1-x0))/d01 * 
                                          std::log((r0+r1+d01) / (r0+r1-d01))
        + ((x-x1)*(y2-y1) - (y-y1)*(x2-x1))/d12 * 
                                          std::log((r1+r2+d12) / (r1+r2-d12))
        + ((x-x2)*(y0-y2) - (y-y2)*(x0-x2))/d20 * 
                                          std::log((r2+r0+d20) / (r2+r0-d20)) );
    
  if (not onpanel)
  {
    phi += z / (4.*M_PI) * (
          std::atan((m01*e0 - h0) / (z*r0)) - std::atan((m01*e1 - h1) / (z*r1))
        + std::atan((m12*e1 - h1) / (z*r1)) - std::atan((m12*e2 - h2) / (z*r2))
        + std::atan((m20*e2 - h2) / (z*r2)) - std::atan((m20*e0 - h0) / (z*r0)) 
        );
  }

  return phi;
}

/******************************************************************************/
//
// Computes potential at a point due to planar quadrilateral source panel with
// unit strength. Panel endpoints given in clockwise order.
//
/******************************************************************************/
double quad_source_potential ( 
                        const double & x, const double & y, const double & z,
                        const double & x0, const double & y0, const double & x1,
                        const double & y1, const double & x2, const double & y2,
                        const double & x3, const double & y3,
                        bool onpanel, const std::string & side )
{
  double d01, d12, d23, d30;
  double r0, r1, r2, r3;
  double dx01, dx12, dx23, dx30, m01, m12, m23, m30;
  double e0, e1, e2, e3;
  double h0, h1, h2, h3;
  double phi, myz;

  if (onpanel)
    myz = 0.;
  else
    myz = z;

  // Geometric quantities needed

  d01 = std::sqrt(std::pow(x1 - x0, 2.) + std::pow(y1 - y0, 2.));
  d12 = std::sqrt(std::pow(x2 - x1, 2.) + std::pow(y2 - y1, 2.));
  d23 = std::sqrt(std::pow(x3 - x2, 2.) + std::pow(y3 - y2, 2.));
  d30 = std::sqrt(std::pow(x0 - x3, 2.) + std::pow(y0 - y3, 2.));

  r0 = std::sqrt(std::pow(x - x0, 2.) + std::pow(y - y0, 2.) + 
                 std::pow(myz, 2.));
  r1 = std::sqrt(std::pow(x - x1, 2.) + std::pow(y - y1, 2.) + 
                 std::pow(myz, 2.));
  r2 = std::sqrt(std::pow(x - x2, 2.) + std::pow(y - y2, 2.) + 
                 std::pow(myz, 2.));
  r3 = std::sqrt(std::pow(x - x3, 2.) + std::pow(y - y3, 2.) + 
                 std::pow(myz, 2.));

  dx01 = sign(x1-x0)*std::max(std::abs(x1-x0), eps);
  dx12 = sign(x2-x1)*std::max(std::abs(x2-x1), eps);
  dx23 = sign(x3-x2)*std::max(std::abs(x3-x2), eps);
  dx30 = sign(x0-x3)*std::max(std::abs(x0-x3), eps);
  m01 = (y1 - y0)/dx01;
  m12 = (y2 - y1)/dx12;
  m23 = (y3 - y2)/dx23;
  m30 = (y0 - y3)/dx30;

  e0 = std::pow(x - x0, 2.) + std::pow(myz, 2.);
  e1 = std::pow(x - x1, 2.) + std::pow(myz, 2.);
  e2 = std::pow(x - x2, 2.) + std::pow(myz, 2.);
  e3 = std::pow(x - x3, 2.) + std::pow(myz, 2.);

  h0 = (x - x0)*(y - y0);
  h1 = (x - x1)*(y - y1);
  h2 = (x - x2)*(y - y2);
  h3 = (x - x3)*(y - y3);

  // Source potential

  phi = -1. / (4.*M_PI) * (
          ((x-x0)*(y1-y0) - (y-y0)*(x1-x0))/d01 * 
                                          std::log((r0+r1+d01) / (r0+r1-d01))
        + ((x-x1)*(y2-y1) - (y-y1)*(x2-x1))/d12 * 
                                          std::log((r1+r2+d12) / (r1+r2-d12))
        + ((x-x2)*(y3-y2) - (y-y2)*(x3-x2))/d23 * 
                                          std::log((r2+r3+d23) / (r2+r3-d23))
        + ((x-x3)*(y0-y3) - (y-y3)*(x0-x3))/d30 * 
                                          std::log((r3+r0+d30) / (r3+r0-d30)) );
    
  if (not onpanel)
  {
    phi += z / (4.*M_PI) * (
          std::atan((m01*e0 - h0) / (z*r0)) - std::atan((m01*e1 - h1) / (z*r1))
        + std::atan((m12*e1 - h1) / (z*r1)) - std::atan((m12*e2 - h2) / (z*r2))
        + std::atan((m23*e2 - h2) / (z*r2)) - std::atan((m23*e3 - h3) / (z*r3))
        + std::atan((m30*e3 - h3) / (z*r3)) - std::atan((m30*e0 - h0) / (z*r0)) 
        );
  }

  return phi;
}

/******************************************************************************/
//
// Computes velocity at a point due to a point source with unit strength. 
// Source is assumed to be at the origin.
//
/******************************************************************************/
Eigen::Vector3d point_source_velocity ( const double & x, const double & y,
                                        const double & z )
{
  double den;
  Eigen::Vector3d velpf;

  den = std::pow(x*x + y*y + z*z, 1.5);

  velpf(0) = x / (4.*M_PI*den);
  velpf(1) = y / (4.*M_PI*den);
  velpf(2) = z / (4.*M_PI*den);

  return velpf;
}

/******************************************************************************/
//
// Computes velocity at a point due to triangular source panel with
// unit strength. Panel endpoints given in clockwise order.
//
/******************************************************************************/
Eigen::Vector3d tri_source_velocity ( 
                        const double & x, const double & y, const double & z,
                        const double & x0, const double & y0, const double & x1,
                        const double & y1, const double & x2, const double & y2,
                        bool onpanel, const std::string & side )
{
  double d01, d12, d20;
  double r0, r1, r2;
  double dx01, dx12, dx20, m01, m12, m20;
  double e0, e1, e2;
  double h0, h1, h2;
  double myz;
  Eigen::Vector3d velpf;

  if (onpanel)
    myz = 0.;
  else
    myz = z;

  // Geometric quantities needed

  d01 = std::sqrt(std::pow(x1 - x0, 2.) + std::pow(y1 - y0, 2.));
  d12 = std::sqrt(std::pow(x2 - x1, 2.) + std::pow(y2 - y1, 2.));
  d20 = std::sqrt(std::pow(x0 - x2, 2.) + std::pow(y0 - y2, 2.));

  r0 = std::sqrt(std::pow(x - x0, 2.) + std::pow(y - y0, 2.) + 
                 std::pow(myz, 2.));
  r1 = std::sqrt(std::pow(x - x1, 2.) + std::pow(y - y1, 2.) + 
                 std::pow(myz, 2.));
  r2 = std::sqrt(std::pow(x - x2, 2.) + std::pow(y - y2, 2.) + 
                 std::pow(myz, 2.));

  if (not onpanel)
  {
    dx01 = sign(x1-x0)*std::max(std::abs(x1-x0), eps);
    dx12 = sign(x2-x1)*std::max(std::abs(x2-x1), eps);
    dx20 = sign(x0-x2)*std::max(std::abs(x0-x2), eps);
    m01 = (y1 - y0)/dx01;
    m12 = (y2 - y1)/dx12;
    m20 = (y0 - y2)/dx20;

    e0 = std::pow(x - x0, 2.) + std::pow(z, 2.);
    e1 = std::pow(x - x1, 2.) + std::pow(z, 2.);
    e2 = std::pow(x - x2, 2.) + std::pow(z, 2.);

    h0 = (x - x0)*(y - y0);
    h1 = (x - x1)*(y - y1);
    h2 = (x - x2)*(y - y2);
  }

  // Velocity components in panel frame

  velpf(0) = 1. / (4.*M_PI) * (
             (y1 - y0)/d01 * std::log((r0 + r1 - d01)/(r0 + r1 + d01))
           + (y2 - y1)/d12 * std::log((r1 + r2 - d12)/(r1 + r2 + d12))
           + (y0 - y2)/d20 * std::log((r2 + r0 - d20)/(r2 + r0 + d20)) );

  velpf(1) = 1. / (4.*M_PI) * (
             (x0 - x1)/d01 * std::log((r0 + r1 - d01)/(r0 + r1 + d01))
           + (x1 - x2)/d12 * std::log((r1 + r2 - d12)/(r1 + r2 + d12))
           + (x2 - x0)/d20 * std::log((r2 + r0 - d20)/(r2 + r0 + d20)) );

  if (not onpanel)
  {
    velpf(2) = 1. / (4.*M_PI) * (
               std::atan((m01*e0 - h0)/(z*r0)) - std::atan((m01*e1 - h1)/(z*r1))
             + std::atan((m12*e1 - h1)/(z*r1)) - std::atan((m12*e2 - h2)/(z*r2))
             + std::atan((m20*e2 - h2)/(z*r2)) - std::atan((m20*e0 - h0)/(z*r0))
             );
  }
  else 
  { 
    if (side == "top") { velpf(2) = 1./2.; }
    else { velpf(2) = -1./2.; }
  }

  return velpf;
}

/******************************************************************************/
//
// Computes velocity at a point due to planar quadrilateral source panel with
// unit strength. Panel endpoints given in clockwise order.
//
/******************************************************************************/
Eigen::Vector3d quad_source_velocity ( 
                        const double & x, const double & y, const double & z,
                        const double & x0, const double & y0, const double & x1,
                        const double & y1, const double & x2, const double & y2,
                        const double & x3, const double & y3, 
                        bool onpanel, const std::string & side )
{
  double d01, d12, d23, d30;
  double r0, r1, r2, r3;
  double dx01, dx12, dx23, dx30, m01, m12, m23, m30;
  double e0, e1, e2, e3;
  double h0, h1, h2, h3;
  double myz;
  Eigen::Vector3d velpf;

  if (onpanel)
    myz = 0.;
  else
    myz = z;

  // Geometric quantities needed

  d01 = std::sqrt(std::pow(x1 - x0, 2.) + std::pow(y1 - y0, 2.));
  d12 = std::sqrt(std::pow(x2 - x1, 2.) + std::pow(y2 - y1, 2.));
  d23 = std::sqrt(std::pow(x3 - x2, 2.) + std::pow(y3 - y2, 2.));
  d30 = std::sqrt(std::pow(x0 - x3, 2.) + std::pow(y0 - y3, 2.));

  r0 = std::sqrt(std::pow(x - x0, 2.) + std::pow(y - y0, 2.) + 
                 std::pow(myz, 2.));
  r1 = std::sqrt(std::pow(x - x1, 2.) + std::pow(y - y1, 2.) + 
                 std::pow(myz, 2.));
  r2 = std::sqrt(std::pow(x - x2, 2.) + std::pow(y - y2, 2.) + 
                 std::pow(myz, 2.));
  r3 = std::sqrt(std::pow(x - x3, 2.) + std::pow(y - y3, 2.) + 
                 std::pow(myz, 2.));

  if (not onpanel)
  {
    dx01 = sign(x1-x0)*std::max(std::abs(x1-x0), eps);
    dx12 = sign(x2-x1)*std::max(std::abs(x2-x1), eps);
    dx23 = sign(x3-x2)*std::max(std::abs(x3-x2), eps);
    dx30 = sign(x0-x3)*std::max(std::abs(x0-x3), eps);
    m01 = (y1 - y0)/dx01;
    m12 = (y2 - y1)/dx12;
    m23 = (y3 - y2)/dx23;
    m30 = (y0 - y3)/dx30;
  
    e0 = std::pow(x - x0, 2.) + std::pow(z, 2.);
    e1 = std::pow(x - x1, 2.) + std::pow(z, 2.);
    e2 = std::pow(x - x2, 2.) + std::pow(z, 2.);
    e3 = std::pow(x - x3, 2.) + std::pow(z, 2.);
  
    h0 = (x - x0)*(y - y0);
    h1 = (x - x1)*(y - y1);
    h2 = (x - x2)*(y - y2);
    h3 = (x - x3)*(y - y3);
  }

  // Velocity components in panel frame

  velpf(0) = 1. / (4.*M_PI) * (
             (y1 - y0)/d01 * std::log((r0 + r1 - d01)/(r0 + r1 + d01))
           + (y2 - y1)/d12 * std::log((r1 + r2 - d12)/(r1 + r2 + d12))
           + (y3 - y2)/d23 * std::log((r2 + r3 - d23)/(r2 + r3 + d23))
           + (y0 - y3)/d30 * std::log((r3 + r0 - d30)/(r3 + r0 + d30)) );

  velpf(1) = 1. / (4.*M_PI) * (
             (x0 - x1)/d01 * std::log((r0 + r1 - d01)/(r0 + r1 + d01))
           + (x1 - x2)/d12 * std::log((r1 + r2 - d12)/(r1 + r2 + d12))
           + (x2 - x3)/d23 * std::log((r2 + r3 - d23)/(r2 + r3 + d23))
           + (x3 - x0)/d30 * std::log((r3 + r0 - d30)/(r3 + r0 + d30)) );

  if (not onpanel)
  {
    velpf(2) = 1. / (4.*M_PI) * (
               std::atan((m01*e0 - h0)/(z*r0)) - std::atan((m01*e1 - h1)/(z*r1))
             + std::atan((m12*e1 - h1)/(z*r1)) - std::atan((m12*e2 - h2)/(z*r2))
             + std::atan((m23*e2 - h2)/(z*r2)) - std::atan((m23*e3 - h3)/(z*r3))
             + std::atan((m30*e3 - h3)/(z*r3)) - std::atan((m30*e0 - h0)/(z*r0))
             );
  }
  else 
  { 
    if (side == "top") { velpf(2) = 1./2.; }
    else { velpf(2) = -1./2.; }
  }

  return velpf;
}

/******************************************************************************/
//
// Computes potential at a point due to a point doublet with unit strength. 
// Doublet is assumed to be at the origin.
//
/******************************************************************************/
double point_doublet_potential ( const double & x, const double & y,
                                 const double & z )
{
  return z / (4.*M_PI*std::pow(x*x + y*y + z*z, 1.5));
}

/******************************************************************************/
//
// Computes potential at a point due to triangular doublet panel with
// unit strength. Panel endpoints given in clockwise order.
//
/******************************************************************************/
double tri_doublet_potential ( 
                        const double & x, const double & y, const double & z,
                        const double & x0, const double & y0, const double & x1,
                        const double & y1, const double & x2, const double & y2,
                        bool onpanel, const std::string & side )
{
  double r0, r1, r2;
  double dx01, dx12, dx20, m01, m12, m20;
  double e0, e1, e2;
  double h0, h1, h2;
  double phi;

  if (onpanel)
  {
    if (side == "top") { phi = -0.5; }
    else { phi = 0.5; }
  }

  else
  {
    // Geometric quantities needed
  
    r0 = std::sqrt(std::pow(x - x0, 2.) + std::pow(y - y0, 2.) + 
                   std::pow(z, 2.));
    r1 = std::sqrt(std::pow(x - x1, 2.) + std::pow(y - y1, 2.) + 
                   std::pow(z, 2.));
    r2 = std::sqrt(std::pow(x - x2, 2.) + std::pow(y - y2, 2.) + 
                   std::pow(z, 2.));
  
    dx01 = sign(x1-x0)*std::max(std::abs(x1-x0), eps);
    dx12 = sign(x2-x1)*std::max(std::abs(x2-x1), eps);
    dx20 = sign(x0-x2)*std::max(std::abs(x0-x2), eps);
    m01 = (y1 - y0)/dx01;
    m12 = (y2 - y1)/dx12;
    m20 = (y0 - y2)/dx20;
  
    e0 = std::pow(x - x0, 2.) + std::pow(z, 2.);
    e1 = std::pow(x - x1, 2.) + std::pow(z, 2.);
    e2 = std::pow(x - x2, 2.) + std::pow(z, 2.);
  
    h0 = (x - x0)*(y - y0);
    h1 = (x - x1)*(y - y1);
    h2 = (x - x2)*(y - y2);
  
    // Doublet potential
  
    phi = 1. / (4.*M_PI) * (
          std::atan((m01*e0 - h0)/(z*r0)) - std::atan((m01*e1 - h1)/(z*r1))
        + std::atan((m12*e1 - h1)/(z*r1)) - std::atan((m12*e2 - h2)/(z*r2))
        + std::atan((m20*e2 - h2)/(z*r2)) - std::atan((m20*e0 - h0)/(z*r0)) );
  }

  return phi;
}

/******************************************************************************/
//
// Computes potential at a point due to planar quadrilateral doublet panel with
// unit strength. Panel endpoints given in clockwise order.
//
/******************************************************************************/
double quad_doublet_potential ( 
                        const double & x, const double & y, const double & z,
                        const double & x0, const double & y0, const double & x1,
                        const double & y1, const double & x2, const double & y2,
                        const double & x3, const double & y3, 
                        bool onpanel, const std::string & side )
{
  double r0, r1, r2, r3;
  double dx01, dx12, dx23, dx30, m01, m12, m23, m30;
  double e0, e1, e2, e3;
  double h0, h1, h2, h3;
  double phi;

  if (onpanel)
  {
    if (side == "top") { phi = -1./2.; }
    else { phi = 1./2.; }
  }

  else
  {
    // Geometric quantities needed
  
    r0 = std::sqrt(std::pow(x - x0, 2.) + std::pow(y - y0, 2.) + 
                   std::pow(z, 2.));
    r1 = std::sqrt(std::pow(x - x1, 2.) + std::pow(y - y1, 2.) + 
                   std::pow(z, 2.));
    r2 = std::sqrt(std::pow(x - x2, 2.) + std::pow(y - y2, 2.) + 
                   std::pow(z, 2.));
    r3 = std::sqrt(std::pow(x - x3, 2.) + std::pow(y - y3, 2.) + 
                   std::pow(z, 2.));
  
    dx01 = sign(x1-x0)*std::max(std::abs(x1-x0), eps);
    dx12 = sign(x2-x1)*std::max(std::abs(x2-x1), eps);
    dx23 = sign(x3-x2)*std::max(std::abs(x3-x2), eps);
    dx30 = sign(x0-x3)*std::max(std::abs(x0-x3), eps);
    m01 = (y1 - y0)/dx01;
    m12 = (y2 - y1)/dx12;
    m23 = (y3 - y2)/dx23;
    m30 = (y0 - y3)/dx30;
  
    e0 = std::pow(x - x0, 2.) + std::pow(z, 2.);
    e1 = std::pow(x - x1, 2.) + std::pow(z, 2.);
    e2 = std::pow(x - x2, 2.) + std::pow(z, 2.);
    e3 = std::pow(x - x3, 2.) + std::pow(z, 2.);
  
    h0 = (x - x0)*(y - y0);
    h1 = (x - x1)*(y - y1);
    h2 = (x - x2)*(y - y2);
    h3 = (x - x3)*(y - y3);
  
    // Doublet potential
  
    phi = 1. / (4.*M_PI) * (
          std::atan((m01*e0 - h0)/(z*r0)) - std::atan((m01*e1 - h1)/(z*r1))
        + std::atan((m12*e1 - h1)/(z*r1)) - std::atan((m12*e2 - h2)/(z*r2))
        + std::atan((m23*e2 - h2)/(z*r2)) - std::atan((m23*e3 - h3)/(z*r3))
        + std::atan((m30*e3 - h3)/(z*r3)) - std::atan((m30*e0 - h0)/(z*r0)) );
  }

  return phi;
}

/******************************************************************************/
//
// Computes velocity at a point due to a point doublet with unit strength. 
// Doublet is assumed to be at the origin.
//
/******************************************************************************/
Eigen::Vector3d point_doublet_velocity ( const double & x, const double & y,
                                         const double & z )
{
  double den;
  Eigen::Vector3d velpf;

  den = std::pow(x*x + y*y + z*z, 2.5);

  velpf(0) = -3.*x*z / (4.*M_PI*den);
  velpf(1) = -3.*y*z / (4.*M_PI*den);
  velpf(2) = (x*x + y*y - 2.*z*z) / (4.*M_PI*den);

  return velpf;
}

/******************************************************************************/
//
// Computes velocity at a point due to triangular doublet panel with
// unit strength. Panel endpoints given in clockwise order.
//
/******************************************************************************/
Eigen::Vector3d tri_doublet_velocity ( 
                        const double & x, const double & y, const double & z,
                        const double & x0, const double & y0, const double & x1,
                        const double & y1, const double & x2, const double & y2,
                        bool onpanel, const std::string & side )
{
  double r0, r1, r2;
  double dx01, dx12, dx20, m01, m12, m20;
  double e0, e1, e2;
  double h0, h1, h2;
  double dr0dx, dr1dx, dr2dx;
  double dr0dy, dr1dy, dr2dy;
  double dr0dz, dr1dz, dr2dz;
  double de0dx, de1dx, de2dx;
  double de0dz, de1dz, de2dz;
  double dh0dx, dh1dx, dh2dx;
  double dh0dy, dh1dy, dh2dy;
  double den0, den1, den2, den3, den4, den5;
  double term1, term2, term3, term4, term5, term6;
  double myz;
  Eigen::Vector3d velpf;

  if (onpanel)
    myz = 0.;
  else
    myz = z;

  // Geometric quantities needed

  r0 = std::sqrt(std::pow(x - x0, 2.) + std::pow(y - y0, 2.) + 
                 std::pow(myz, 2.));
  r1 = std::sqrt(std::pow(x - x1, 2.) + std::pow(y - y1, 2.) + 
                 std::pow(myz, 2.));
  r2 = std::sqrt(std::pow(x - x2, 2.) + std::pow(y - y2, 2.) + 
                 std::pow(myz, 2.));

  dx01 = sign(x1-x0)*std::max(std::abs(x1-x0), eps);
  dx12 = sign(x2-x1)*std::max(std::abs(x2-x1), eps);
  dx20 = sign(x0-x2)*std::max(std::abs(x0-x2), eps);
  m01 = (y1 - y0)/dx01;
  m12 = (y2 - y1)/dx12;
  m20 = (y0 - y2)/dx20;

  e0 = std::pow(x - x0, 2.) + std::pow(myz, 2.);
  e1 = std::pow(x - x1, 2.) + std::pow(myz, 2.);
  e2 = std::pow(x - x2, 2.) + std::pow(myz, 2.);

  h0 = (x - x0)*(y - y0);
  h1 = (x - x1)*(y - y1);
  h2 = (x - x2)*(y - y2);

  if (not onpanel)
  {
    // Derivatives of geometric quantities
  
    dr0dx = (x - x0)/r0;
    dr1dx = (x - x1)/r1;
    dr2dx = (x - x2)/r2;
  
    dr0dy = (y - y0)/r0;
    dr1dy = (y - y1)/r1;
    dr2dy = (y - y2)/r2;
  
    dr0dz = z/r0;
    dr1dz = z/r1;
    dr2dz = z/r2;
  
    de0dx = 2.*(x - x0);
    de1dx = 2.*(x - x1);
    de2dx = 2.*(x - x2);
  
    de0dz = 2.*z;
    de1dz = 2.*z;
    de2dz = 2.*z;
  
    dh0dx = (y - y0);
    dh1dx = (y - y1);
    dh2dx = (y - y2);
  
    dh0dy = (x - x0);
    dh1dy = (x - x1);
    dh2dy = (x - x2);
  
    // Denominator terms
  
    den0 = std::pow(z*r0, 2.) + std::pow(m01*e0 - h0, 2.);
    den1 = std::pow(z*r1, 2.) + std::pow(m01*e1 - h1, 2.);
    den2 = std::pow(z*r1, 2.) + std::pow(m12*e1 - h1, 2.);
    den3 = std::pow(z*r2, 2.) + std::pow(m12*e2 - h2, 2.);
    den4 = std::pow(z*r2, 2.) + std::pow(m20*e2 - h2, 2.);
    den5 = std::pow(z*r0, 2.) + std::pow(m20*e0 - h0, 2.);
  
    // Velocity components in panel frame
  
    velpf(0) = 1. / (4.*M_PI) * z * (
               (r0*(m01*de0dx - dh0dx) - dr0dx*(m01*e0 - h0)) / den0
             - (r1*(m01*de1dx - dh1dx) - dr1dx*(m01*e1 - h1)) / den1
             + (r1*(m12*de1dx - dh1dx) - dr1dx*(m12*e1 - h1)) / den2
             - (r2*(m12*de2dx - dh2dx) - dr2dx*(m12*e2 - h2)) / den3
             + (r2*(m20*de2dx - dh2dx) - dr2dx*(m20*e2 - h2)) / den4
             - (r0*(m20*de0dx - dh0dx) - dr0dx*(m20*e0 - h0)) / den5 );
  
    velpf(1) = 1. / (4.*M_PI) * z * (
               (-r0*dh0dy - dr0dy*(m01*e0 - h0)) / den0
             - (-r1*dh1dy - dr1dy*(m01*e1 - h1)) / den1
             + (-r1*dh1dy - dr1dy*(m12*e1 - h1)) / den2
             - (-r2*dh2dy - dr2dy*(m12*e2 - h2)) / den3
             + (-r2*dh2dy - dr2dy*(m20*e2 - h2)) / den4
             - (-r0*dh0dy - dr0dy*(m20*e0 - h0)) / den5 );
  
    velpf(2) = 1. / (4.*M_PI) * (
               (z*r0*m01*de0dz - (z*dr0dz + r0)*(m01*e0 - h0)) / den0
             - (z*r1*m01*de1dz - (z*dr1dz + r1)*(m01*e1 - h1)) / den1
             + (z*r1*m12*de1dz - (z*dr1dz + r1)*(m12*e1 - h1)) / den2
             - (z*r2*m12*de2dz - (z*dr2dz + r2)*(m12*e2 - h2)) / den3
             + (z*r2*m20*de2dz - (z*dr2dz + r2)*(m20*e2 - h2)) / den4
             - (z*r0*m20*de0dz - (z*dr0dz + r0)*(m20*e0 - h0)) / den5 );
  }
  
  else
  {
    velpf(0) = 0.0;
    velpf(1) = 0.0;

    term1 = m01*e0 - h0;
    term2 = m01*e1 - h1;
    term3 = m12*e1 - h1;
    term4 = m12*e2 - h2;
    term5 = m20*e2 - h2;
    term6 = m20*e0 - h0;
    velpf(2) = 1. / (4.*M_PI) * (
               -r0/term1 + r1/term2 
             -  r1/term3 + r2/term4
             -  r2/term5 + r0/term6 );
  }

  return velpf;
}

/******************************************************************************/
//
// Computes velocity at a point due to planar quadrilateral doublet panel with
// unit strength. Panel endpoints given in clockwise order.
//
/******************************************************************************/
Eigen::Vector3d quad_doublet_velocity ( 
                        const double & x, const double & y, const double & z,
                        const double & x0, const double & y0, const double & x1,
                        const double & y1, const double & x2, const double & y2,
                        const double & x3, const double & y3, 
                        bool onpanel, const std::string & side )
{
  double r0, r1, r2, r3;
  double dx01, dx12, dx23, dx30, m01, m12, m23, m30;
  double e0, e1, e2, e3;
  double h0, h1, h2, h3;
  double dr0dx, dr1dx, dr2dx, dr3dx;
  double dr0dy, dr1dy, dr2dy, dr3dy;
  double dr0dz, dr1dz, dr2dz, dr3dz;
  double de0dx, de1dx, de2dx, de3dx;
  double de0dz, de1dz, de2dz, de3dz;
  double dh0dx, dh1dx, dh2dx, dh3dx;
  double dh0dy, dh1dy, dh2dy, dh3dy;
  double den0, den1, den2, den3, den4, den5, den6, den7;
  double term1, term2, term3, term4, term5, term6, term7, term8;
  double myz;
  Eigen::Vector3d velpf;

  if (onpanel)
    myz = 0.;
  else
    myz = z;

  // Geometric quantities needed

  r0 = std::sqrt(std::pow(x - x0, 2.) + std::pow(y - y0, 2.) + 
                 std::pow(myz, 2.));
  r1 = std::sqrt(std::pow(x - x1, 2.) + std::pow(y - y1, 2.) + 
                 std::pow(myz, 2.));
  r2 = std::sqrt(std::pow(x - x2, 2.) + std::pow(y - y2, 2.) + 
                 std::pow(myz, 2.));
  r3 = std::sqrt(std::pow(x - x3, 2.) + std::pow(y - y3, 2.) + 
                 std::pow(myz, 2.));

  dx01 = sign(x1-x0)*std::max(std::abs(x1-x0), eps);
  dx12 = sign(x2-x1)*std::max(std::abs(x2-x1), eps);
  dx23 = sign(x3-x2)*std::max(std::abs(x3-x2), eps);
  dx30 = sign(x0-x3)*std::max(std::abs(x0-x3), eps);
  m01 = (y1 - y0)/dx01;
  m12 = (y2 - y1)/dx12;
  m23 = (y3 - y2)/dx23;
  m30 = (y0 - y3)/dx30;

  e0 = std::pow(x - x0, 2.) + std::pow(myz, 2.);
  e1 = std::pow(x - x1, 2.) + std::pow(myz, 2.);
  e2 = std::pow(x - x2, 2.) + std::pow(myz, 2.);
  e3 = std::pow(x - x3, 2.) + std::pow(myz, 2.);

  h0 = (x - x0)*(y - y0);
  h1 = (x - x1)*(y - y1);
  h2 = (x - x2)*(y - y2);
  h3 = (x - x3)*(y - y3);

  if (not onpanel)
  {
    // Derivatives of geometric quantities
  
    dr0dx = (x - x0)/r0;
    dr1dx = (x - x1)/r1;
    dr2dx = (x - x2)/r2;
    dr3dx = (x - x3)/r3;
  
    dr0dy = (y - y0)/r0;
    dr1dy = (y - y1)/r1;
    dr2dy = (y - y2)/r2;
    dr3dy = (y - y3)/r3;
  
    dr0dz = z/r0;
    dr1dz = z/r1;
    dr2dz = z/r2;
    dr3dz = z/r3;
  
    de0dx = 2.*(x - x0);
    de1dx = 2.*(x - x1);
    de2dx = 2.*(x - x2);
    de3dx = 2.*(x - x3);
  
    de0dz = 2.*z;
    de1dz = 2.*z;
    de2dz = 2.*z;
    de3dz = 2.*z;
  
    dh0dx = (y - y0);
    dh1dx = (y - y1);
    dh2dx = (y - y2);
    dh3dx = (y - y3);
  
    dh0dy = (x - x0);
    dh1dy = (x - x1);
    dh2dy = (x - x2);
    dh3dy = (x - x3);
  
    // Denominator terms
  
    den0 = std::pow(z*r0, 2.) + std::pow(m01*e0 - h0, 2.);
    den1 = std::pow(z*r1, 2.) + std::pow(m01*e1 - h1, 2.);
    den2 = std::pow(z*r1, 2.) + std::pow(m12*e1 - h1, 2.);
    den3 = std::pow(z*r2, 2.) + std::pow(m12*e2 - h2, 2.);
    den4 = std::pow(z*r2, 2.) + std::pow(m23*e2 - h2, 2.);
    den5 = std::pow(z*r3, 2.) + std::pow(m23*e3 - h3, 2.);
    den6 = std::pow(z*r3, 2.) + std::pow(m30*e3 - h3, 2.);
    den7 = std::pow(z*r0, 2.) + std::pow(m30*e0 - h0, 2.);
  
    // Velocity components in panel frame
  
    velpf(0) = 1. / (4.*M_PI) * z * (
               (r0*(m01*de0dx - dh0dx) - dr0dx*(m01*e0 - h0)) / den0
             - (r1*(m01*de1dx - dh1dx) - dr1dx*(m01*e1 - h1)) / den1
             + (r1*(m12*de1dx - dh1dx) - dr1dx*(m12*e1 - h1)) / den2
             - (r2*(m12*de2dx - dh2dx) - dr2dx*(m12*e2 - h2)) / den3
             + (r2*(m23*de2dx - dh2dx) - dr2dx*(m23*e2 - h2)) / den4
             - (r3*(m23*de3dx - dh3dx) - dr3dx*(m23*e3 - h3)) / den5
             + (r3*(m30*de3dx - dh3dx) - dr3dx*(m30*e3 - h3)) / den6
             - (r0*(m30*de0dx - dh0dx) - dr0dx*(m30*e0 - h0)) / den7 );
  
    velpf(1) = 1. / (4.*M_PI) * z * (
               (-r0*dh0dy - dr0dy*(m01*e0 - h0)) / den0
             - (-r1*dh1dy - dr1dy*(m01*e1 - h1)) / den1
             + (-r1*dh1dy - dr1dy*(m12*e1 - h1)) / den2
             - (-r2*dh2dy - dr2dy*(m12*e2 - h2)) / den3
             + (-r2*dh2dy - dr2dy*(m23*e2 - h2)) / den4
             - (-r3*dh3dy - dr3dy*(m23*e3 - h3)) / den5
             + (-r3*dh3dy - dr3dy*(m30*e3 - h3)) / den6
             - (-r0*dh0dy - dr0dy*(m30*e0 - h0)) / den7 );
  
    velpf(2) = 1. / (4.*M_PI) * (
               (z*r0*m01*de0dz - (z*dr0dz + r0)*(m01*e0 - h0)) / den0
             - (z*r1*m01*de1dz - (z*dr1dz + r1)*(m01*e1 - h1)) / den1
             + (z*r1*m12*de1dz - (z*dr1dz + r1)*(m12*e1 - h1)) / den2
             - (z*r2*m12*de2dz - (z*dr2dz + r2)*(m12*e2 - h2)) / den3
             + (z*r2*m23*de2dz - (z*dr2dz + r2)*(m23*e2 - h2)) / den4
             - (z*r3*m23*de3dz - (z*dr3dz + r3)*(m23*e3 - h3)) / den5
             + (z*r3*m30*de3dz - (z*dr3dz + r3)*(m30*e3 - h3)) / den6
             - (z*r0*m30*de0dz - (z*dr0dz + r0)*(m30*e0 - h0)) / den7 );
  }
  
  else
  {
    velpf(0) = 0.0;
    velpf(1) = 0.0;

    term1 = m01*e0 - h0;
    term2 = m01*e1 - h1;
    term3 = m12*e1 - h1;
    term4 = m12*e2 - h2;
    term5 = m23*e2 - h2; 
    term6 = m23*e3 - h3;
    term7 = m30*e3 - h3;
    term8 = m30*e0 - h0;
    velpf(2) = 1. / (4.*M_PI) * (
               -r0/term1 + r1/term2 
             -  r1/term3 + r2/term4 
             -  r2/term5 + r3/term6
             -  r3/term7 + r0/term8 );
  }

  return velpf;
}

/******************************************************************************/
//
// Computes velocity at a point due to a vortex filament of finite or semi-
// infinite length with finite core radius and unit strength
//
/******************************************************************************/
Eigen::Vector3d vortex_velocity ( 
                        const double & x, const double & y, const double & z,
                        const double & x1, const double & y1, const double & z1,
                        const double & x2, const double & y2, const double & z2,
                        const double & rcore, bool semi_infinite )
{
  double l0, l1, l2, l1x2, d, f1, cosB1, cosB2, r0dr1, r0dr2, term;
  Eigen::Vector3d r0, r1, r2, r1xr2, vel;

  // Vector from point 1 to point 2 (and length)

  r0(0) = x2 - x1;
  r0(1) = y2 - y1;
  r0(2) = z2 - z1;
  l0 = r0.norm();

  if (std::abs(l0) < eps)
  {
    vel(0) = 0.0;
    vel(1) = 0.0;
    vel(2) = 0.0;
    return vel;
  }  

  // Vectors from endpoints to point P (and lengths)

  r1(0) = x - x1;
  r1(1) = y - y1;
  r1(2) = z - z1;
  l1 = r1.norm();

  r2(0) = x - x2;
  r2(1) = y - y2;
  r2(2) = z - z2;
  l2 = r2.norm();

  // Cross product of r1 and r2 (and length)

  r1xr2(0) = r1(1)*r2(2) - r1(2)*r2(1);
  r1xr2(1) = r1(2)*r2(0) - r1(0)*r2(2);
  r1xr2(2) = r1(0)*r2(1) - r1(1)*r2(0);
  l1x2 = r1xr2.norm();

  // Distance to point P

  d = l1x2 / l0;
  if (d < eps)
  {
    vel(0) = 0.0;
    vel(1) = 0.0;
    vel(2) = 0.0;
    return vel;
  }  

  // Finite core model approximating 1/d outside of rcore

  f1 = d / std::pow(d + corefact*rcore, 2.0);

  // Dot products r0*r1 and r0*r2

  r0dr1 = r0(0)*r1(0) + r0(1)*r1(1) + r0(2)*r1(2);
  r0dr2 = r0(0)*r2(0) + r0(1)*r2(1) + r0(2)*r2(2);
  
  // Cosines of angles

  cosB1 = r0dr1 / (l0*l1);
  if (semi_infinite) { cosB2 = -1.; }
  else { cosB2 = r0dr2 / (l0*l2); }

  // Induced velocity

  term = f1 / (4.*M_PI*l1x2) * (cosB1 - cosB2);
  vel(0) = term * r1xr2(0);
  vel(1) = term * r1xr2(1);
  vel(2) = term * r1xr2(2);

  return vel;
}
